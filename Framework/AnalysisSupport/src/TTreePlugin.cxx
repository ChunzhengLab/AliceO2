// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Framework/RootArrowFilesystem.h"
#include "Framework/Plugins.h"
#include "Framework/Signpost.h"
#include "Framework/Endian.h"
#include <arrow/dataset/file_base.h>
#include <arrow/util/key_value_metadata.h>
#include <arrow/array/array_nested.h>
#include <arrow/array/array_primitive.h>
#include <arrow/array/builder_nested.h>
#include <arrow/array/builder_primitive.h>
#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>
#include <TLeaf.h>

O2_DECLARE_DYNAMIC_LOG(root_arrow_fs);

namespace o2::framework
{

class TTreeFileWriteOptions : public arrow::dataset::FileWriteOptions
{
 public:
  TTreeFileWriteOptions(std::shared_ptr<arrow::dataset::FileFormat> format)
    : FileWriteOptions(format)
  {
  }
};

// A filesystem which allows me to get a TTree
class TTreeFileSystem : public VirtualRootFileSystemBase
{
 public:
  ~TTreeFileSystem() override;

  std::shared_ptr<VirtualRootFileSystemBase> GetSubFilesystem(arrow::dataset::FileSource source) override
  {
    return std::dynamic_pointer_cast<VirtualRootFileSystemBase>(shared_from_this());
  };

  arrow::Result<std::shared_ptr<arrow::io::OutputStream>> OpenOutputStream(
    const std::string& path,
    const std::shared_ptr<const arrow::KeyValueMetadata>& metadata) override;

  virtual std::unique_ptr<TTree>& GetTree(arrow::dataset::FileSource source) = 0;
};

class SingleTreeFileSystem : public TTreeFileSystem
{
 public:
  SingleTreeFileSystem(TTree* tree)
    : TTreeFileSystem(),
      mTree(tree)
  {
  }

  arrow::Result<arrow::fs::FileInfo> GetFileInfo(std::string const& path) override;

  std::string type_name() const override
  {
    return "ttree";
  }

  std::unique_ptr<TTree>& GetTree(arrow::dataset::FileSource) override
  {
    // Simply return the only TTree we have
    return mTree;
  }

 private:
  std::unique_ptr<TTree> mTree;
};

arrow::Result<arrow::fs::FileInfo> SingleTreeFileSystem::GetFileInfo(std::string const& path)
{
  arrow::dataset::FileSource source(path, shared_from_this());
  arrow::fs::FileInfo result;
  result.set_path(path);
  result.set_type(arrow::fs::FileType::File);
  return result;
}

class TTreeFileFragment : public arrow::dataset::FileFragment
{
 public:
  TTreeFileFragment(arrow::dataset::FileSource source,
                    std::shared_ptr<arrow::dataset::FileFormat> format,
                    arrow::compute::Expression partition_expression,
                    std::shared_ptr<arrow::Schema> physical_schema)
    : FileFragment(std::move(source), std::move(format), std::move(partition_expression), std::move(physical_schema))
  {
  }
};

class TTreeFileFormat : public arrow::dataset::FileFormat
{
  size_t& mTotCompressedSize;
  size_t& mTotUncompressedSize;

 public:
  TTreeFileFormat(size_t& totalCompressedSize, size_t& totalUncompressedSize)
    : FileFormat({}),
      mTotCompressedSize(totalCompressedSize),
      mTotUncompressedSize(totalUncompressedSize)
  {
  }

  ~TTreeFileFormat() override = default;

  std::string type_name() const override
  {
    return "ttree";
  }

  bool Equals(const FileFormat& other) const override
  {
    return other.type_name() == this->type_name();
  }

  arrow::Result<bool> IsSupported(const arrow::dataset::FileSource& source) const override
  {
    auto fs = std::dynamic_pointer_cast<VirtualRootFileSystemBase>(source.filesystem());
    auto subFs = fs->GetSubFilesystem(source);
    if (std::dynamic_pointer_cast<TTreeFileSystem>(subFs)) {
      return true;
    }
    return false;
  }

  arrow::Result<std::shared_ptr<arrow::Schema>> Inspect(const arrow::dataset::FileSource& source) const override;
  /// \brief Create a FileFragment for a FileSource.
  arrow::Result<std::shared_ptr<arrow::dataset::FileFragment>> MakeFragment(
    arrow::dataset::FileSource source, arrow::compute::Expression partition_expression,
    std::shared_ptr<arrow::Schema> physical_schema) override;

  arrow::Result<std::shared_ptr<arrow::dataset::FileWriter>> MakeWriter(std::shared_ptr<arrow::io::OutputStream> destination, std::shared_ptr<arrow::Schema> schema, std::shared_ptr<arrow::dataset::FileWriteOptions> options, arrow::fs::FileLocator destination_locator) const override;

  std::shared_ptr<arrow::dataset::FileWriteOptions> DefaultWriteOptions() override;

  arrow::Result<arrow::RecordBatchGenerator> ScanBatchesAsync(
    const std::shared_ptr<arrow::dataset::ScanOptions>& options,
    const std::shared_ptr<arrow::dataset::FileFragment>& fragment) const override;
};

// An arrow outputstream which allows to write to a TTree. Eventually
// with a prefix for the branches.
class TTreeOutputStream : public arrow::io::OutputStream
{
 public:
  // Using a pointer means that the tree itself is owned by another 
  // class
  TTreeOutputStream(TTree *, std::string branchPrefix);

  arrow::Status Close() override;

  arrow::Result<int64_t> Tell() const override;

  arrow::Status Write(const void* data, int64_t nbytes) override;

  bool closed() const override;

  TBranch* CreateBranch(char const* branchName, char const* sizeBranch);

  TTree* GetTree()
  {
    return mTree;
  }

 private:
  TTree* mTree;
  std::string mBranchPrefix;
};

// An arrow outputstream which allows to write to a ttree
// @a branch prefix is to be used to identify a set of branches which all belong to
// the same table.
TTreeOutputStream::TTreeOutputStream(TTree* f, std::string branchPrefix)
  : mTree(f),
    mBranchPrefix(std::move(branchPrefix))
{
}

arrow::Status TTreeOutputStream::Close()
{
  if (mTree->GetCurrentFile() == nullptr) {
    return arrow::Status::Invalid("Cannot close a tree not attached to a file");
  }
  mTree->GetCurrentFile()->Close();
  return arrow::Status::OK();
}

arrow::Result<int64_t> TTreeOutputStream::Tell() const
{
  return arrow::Result<int64_t>(arrow::Status::NotImplemented("Cannot move"));
}

arrow::Status TTreeOutputStream::Write(const void* data, int64_t nbytes)
{
  return arrow::Status::NotImplemented("Cannot write raw bytes to a TTree");
}

bool TTreeOutputStream::closed() const
{
  // A standalone tree is never closed.
  if (mTree->GetCurrentFile() == nullptr) {
    return false;
  }
  return mTree->GetCurrentFile()->IsOpen() == false;
}

TBranch* TTreeOutputStream::CreateBranch(char const* branchName, char const* sizeBranch)
{
  return mTree->Branch((mBranchPrefix + "/" + branchName).c_str(), (char*)nullptr, (mBranchPrefix + sizeBranch).c_str());
}

struct TTreePluginContext {
  size_t totalCompressedSize = 0;
  size_t totalUncompressedSize = 0;
  std::shared_ptr<o2::framework::TTreeFileFormat> format = nullptr;
};

struct TTreeObjectReadingImplementation : public RootArrowFactoryPlugin {
  RootArrowFactory* create() override
  {
    auto context = new TTreePluginContext;
    context->format = std::make_shared<o2::framework::TTreeFileFormat>(context->totalCompressedSize, context->totalUncompressedSize);
    return new RootArrowFactory{
      .options = [context]() { return context->format->DefaultWriteOptions(); },
      .format = [context]() { return context->format; },
      .getSubFilesystem = [](void* handle) {
        auto tree = (TTree*)handle;
        return std::shared_ptr<VirtualRootFileSystemBase>(new SingleTreeFileSystem(tree)); },
    };
  }
};

arrow::Result<arrow::RecordBatchGenerator> TTreeFileFormat::ScanBatchesAsync(
  const std::shared_ptr<arrow::dataset::ScanOptions>& options,
  const std::shared_ptr<arrow::dataset::FileFragment>& fragment) const
{
  // Get the fragment as a TTreeFragment. This might be PART of a TTree.
  auto treeFragment = std::dynamic_pointer_cast<TTreeFileFragment>(fragment);
  // This is the schema we want to read
  auto dataset_schema = options->dataset_schema;

  auto generator = [pool = options->pool, treeFragment, dataset_schema, &totalCompressedSize = mTotCompressedSize,
                    &totalUncompressedSize = mTotUncompressedSize]() -> arrow::Future<std::shared_ptr<arrow::RecordBatch>> {
    auto schema = treeFragment->format()->Inspect(treeFragment->source());

    std::vector<std::shared_ptr<arrow::Array>> columns;
    std::vector<std::shared_ptr<arrow::Field>> fields = dataset_schema->fields();
    auto physical_schema = *treeFragment->ReadPhysicalSchema();

    static TBufferFile buffer{TBuffer::EMode::kWrite, 4 * 1024 * 1024};
    auto containerFS = std::dynamic_pointer_cast<VirtualRootFileSystemBase>(treeFragment->source().filesystem());
    auto fs = std::dynamic_pointer_cast<TTreeFileSystem>(containerFS->GetSubFilesystem(treeFragment->source()));

    int64_t rows = -1;
    auto& tree = fs->GetTree(treeFragment->source());
    for (auto& field : fields) {
      // The field actually on disk
      auto physicalField = physical_schema->GetFieldByName(field->name());
      TBranch* branch = tree->GetBranch(physicalField->name().c_str());
      assert(branch);
      buffer.Reset();
      auto totalEntries = branch->GetEntries();
      if (rows == -1) {
        rows = totalEntries;
      }
      if (rows != totalEntries) {
        throw runtime_error_f("Unmatching number of rows for branch %s", branch->GetName());
      }
      arrow::Status status;
      int readEntries = 0;
      std::shared_ptr<arrow::Array> array;
      auto listType = std::dynamic_pointer_cast<arrow::FixedSizeListType>(physicalField->type());
      if (physicalField->type() == arrow::boolean() ||
          (listType && physicalField->type()->field(0)->type() == arrow::boolean())) {
        if (listType) {
          std::unique_ptr<arrow::ArrayBuilder> builder = nullptr;
          auto status = arrow::MakeBuilder(pool, physicalField->type()->field(0)->type(), &builder);
          if (!status.ok()) {
            throw runtime_error("Cannot create value builder");
          }
          auto listBuilder = std::make_unique<arrow::FixedSizeListBuilder>(pool, std::move(builder), listType->list_size());
          auto valueBuilder = listBuilder.get()->value_builder();
          // boolean array special case: we need to use builder to create the bitmap
          status = valueBuilder->Reserve(totalEntries * listType->list_size());
          status &= listBuilder->Reserve(totalEntries);
          if (!status.ok()) {
            throw runtime_error("Failed to reserve memory for array builder");
          }
          while (readEntries < totalEntries) {
            auto readLast = branch->GetBulkRead().GetBulkEntries(readEntries, buffer);
            readEntries += readLast;
            status &= static_cast<arrow::BooleanBuilder*>(valueBuilder)->AppendValues(reinterpret_cast<uint8_t const*>(buffer.GetCurrent()), readLast * listType->list_size());
          }
          status &= static_cast<arrow::FixedSizeListBuilder*>(listBuilder.get())->AppendValues(readEntries);
          if (!status.ok()) {
            throw runtime_error("Failed to append values to array");
          }
          status &= listBuilder->Finish(&array);
          if (!status.ok()) {
            throw runtime_error("Failed to create array");
          }
        } else if (listType == nullptr) {
          std::unique_ptr<arrow::ArrayBuilder> builder = nullptr;
          auto status = arrow::MakeBuilder(pool, physicalField->type(), &builder);
          if (!status.ok()) {
            throw runtime_error("Cannot create builder");
          }
          auto valueBuilder = static_cast<arrow::BooleanBuilder*>(builder.get());
          // boolean array special case: we need to use builder to create the bitmap
          status = valueBuilder->Reserve(totalEntries);
          if (!status.ok()) {
            throw runtime_error("Failed to reserve memory for array builder");
          }
          while (readEntries < totalEntries) {
            auto readLast = branch->GetBulkRead().GetBulkEntries(readEntries, buffer);
            readEntries += readLast;
            status &= valueBuilder->AppendValues(reinterpret_cast<uint8_t const*>(buffer.GetCurrent()), readLast);
          }
          if (!status.ok()) {
            throw runtime_error("Failed to append values to array");
          }
          status &= valueBuilder->Finish(&array);
          if (!status.ok()) {
            throw runtime_error("Failed to create array");
          }
        }
      } else {
        // other types: use serialized read to build arrays directly.
        auto typeSize = physicalField->type()->byte_width();
        // This is needed for branches which have not been persisted.
        auto bytes = branch->GetTotBytes();
        auto branchSize = bytes ? bytes : 1000000;
        auto&& result = arrow::AllocateResizableBuffer(branchSize, pool);
        if (!result.ok()) {
          throw runtime_error("Cannot allocate values buffer");
        }
        std::shared_ptr<arrow::Buffer> arrowValuesBuffer = std::move(result).ValueUnsafe();
        auto ptr = arrowValuesBuffer->mutable_data();
        if (ptr == nullptr) {
          throw runtime_error("Invalid buffer");
        }

        std::unique_ptr<TBufferFile> offsetBuffer = nullptr;

        uint32_t offset = 0;
        int count = 0;
        std::shared_ptr<arrow::Buffer> arrowOffsetBuffer;
        std::span<int> offsets;
        int size = 0;
        uint32_t totalSize = 0;
        TBranch* mSizeBranch = nullptr;
        int64_t listSize = 1;
        if (auto fixedSizeList = std::dynamic_pointer_cast<arrow::FixedSizeListType>(physicalField->type())) {
          listSize = fixedSizeList->list_size();
          typeSize = fixedSizeList->field(0)->type()->byte_width();
        } else if (auto vlaListType = std::dynamic_pointer_cast<arrow::ListType>(physicalField->type())) {
          listSize = -1;
          typeSize = vlaListType->field(0)->type()->byte_width();
        }
        if (listSize == -1) {
          mSizeBranch = branch->GetTree()->GetBranch((std::string{branch->GetName()} + "_size").c_str());
          offsetBuffer = std::make_unique<TBufferFile>(TBuffer::EMode::kWrite, 4 * 1024 * 1024);
          result = arrow::AllocateResizableBuffer((totalEntries + 1) * (int64_t)sizeof(int), pool);
          if (!result.ok()) {
            throw runtime_error("Cannot allocate offset buffer");
          }
          arrowOffsetBuffer = std::move(result).ValueUnsafe();
          unsigned char* ptrOffset = arrowOffsetBuffer->mutable_data();
          auto* tPtrOffset = reinterpret_cast<int*>(ptrOffset);
          offsets = std::span<int>{tPtrOffset, tPtrOffset + totalEntries + 1};

          // read sizes first
          while (readEntries < totalEntries) {
            auto readLast = mSizeBranch->GetBulkRead().GetEntriesSerialized(readEntries, *offsetBuffer);
            readEntries += readLast;
            for (auto i = 0; i < readLast; ++i) {
              offsets[count++] = (int)offset;
              offset += swap32_(reinterpret_cast<uint32_t*>(offsetBuffer->GetCurrent())[i]);
            }
          }
          offsets[count] = (int)offset;
          totalSize = offset;
          readEntries = 0;
        }

        while (readEntries < totalEntries) {
          auto readLast = branch->GetBulkRead().GetEntriesSerialized(readEntries, buffer);
          if (listSize == -1) {
            size = offsets[readEntries + readLast] - offsets[readEntries];
          } else {
            size = readLast * listSize;
          }
          readEntries += readLast;
          swapCopy(ptr, buffer.GetCurrent(), size, typeSize);
          ptr += (ptrdiff_t)(size * typeSize);
        }
        if (listSize >= 1) {
          totalSize = readEntries * listSize;
        }
        std::shared_ptr<arrow::PrimitiveArray> varray;
        switch (listSize) {
          case -1:
            varray = std::make_shared<arrow::PrimitiveArray>(physicalField->type()->field(0)->type(), totalSize, arrowValuesBuffer);
            array = std::make_shared<arrow::ListArray>(physicalField->type(), readEntries, arrowOffsetBuffer, varray);
            break;
          case 1:
            array = std::make_shared<arrow::PrimitiveArray>(physicalField->type(), readEntries, arrowValuesBuffer);
            break;
          default:
            varray = std::make_shared<arrow::PrimitiveArray>(physicalField->type()->field(0)->type(), totalSize, arrowValuesBuffer);
            array = std::make_shared<arrow::FixedSizeListArray>(physicalField->type(), readEntries, varray);
        }
      }

      branch->SetStatus(false);
      branch->DropBaskets("all");
      branch->Reset();
      branch->GetTransientBuffer(0)->Expand(0);

      columns.push_back(array);
    }
    auto batch = arrow::RecordBatch::Make(dataset_schema, rows, columns);
    totalCompressedSize += tree->GetZipBytes();
    totalUncompressedSize += tree->GetTotBytes();
    return batch;
  };
  return generator;
}

char const* rootSuffixFromArrow(arrow::Type::type id)
{
  switch (id) {
    case arrow::Type::BOOL:
      return "/O";
    case arrow::Type::UINT8:
      return "/b";
    case arrow::Type::UINT16:
      return "/s";
    case arrow::Type::UINT32:
      return "/i";
    case arrow::Type::UINT64:
      return "/l";
    case arrow::Type::INT8:
      return "/B";
    case arrow::Type::INT16:
      return "/S";
    case arrow::Type::INT32:
      return "/I";
    case arrow::Type::INT64:
      return "/L";
    case arrow::Type::FLOAT:
      return "/F";
    case arrow::Type::DOUBLE:
      return "/D";
    default:
      throw runtime_error("Unsupported arrow column type");
  }
}

arrow::Result<std::shared_ptr<arrow::io::OutputStream>> TTreeFileSystem::OpenOutputStream(
  const std::string& path,
  const std::shared_ptr<const arrow::KeyValueMetadata>& metadata)
{
  arrow::dataset::FileSource source{path, shared_from_this()};
  auto prefix = metadata->Get("branch_prefix");
  if (prefix.ok()) {
    return std::make_shared<TTreeOutputStream>(GetTree(source).get(), *prefix);
  }
  return std::make_shared<TTreeOutputStream>(GetTree(source).get(), "");
}

namespace
{
struct BranchInfo {
  std::string name;
  TBranch* ptr;
  bool mVLA;
};
} // namespace

auto arrowTypeFromROOT(EDataType type, int size)
{
  auto typeGenerator = [](std::shared_ptr<arrow::DataType> const& type, int size) -> std::shared_ptr<arrow::DataType> {
    switch (size) {
      case -1:
        return arrow::list(type);
      case 1:
        return std::move(type);
      default:
        return arrow::fixed_size_list(type, size);
    }
  };

  switch (type) {
    case EDataType::kBool_t:
      return typeGenerator(arrow::boolean(), size);
    case EDataType::kUChar_t:
      return typeGenerator(arrow::uint8(), size);
    case EDataType::kUShort_t:
      return typeGenerator(arrow::uint16(), size);
    case EDataType::kUInt_t:
      return typeGenerator(arrow::uint32(), size);
    case EDataType::kULong64_t:
      return typeGenerator(arrow::uint64(), size);
    case EDataType::kChar_t:
      return typeGenerator(arrow::int8(), size);
    case EDataType::kShort_t:
      return typeGenerator(arrow::int16(), size);
    case EDataType::kInt_t:
      return typeGenerator(arrow::int32(), size);
    case EDataType::kLong64_t:
      return typeGenerator(arrow::int64(), size);
    case EDataType::kFloat_t:
      return typeGenerator(arrow::float32(), size);
    case EDataType::kDouble_t:
      return typeGenerator(arrow::float64(), size);
    default:
      throw o2::framework::runtime_error_f("Unsupported branch type: %d", static_cast<int>(type));
  }
}

arrow::Result<std::shared_ptr<arrow::Schema>> TTreeFileFormat::Inspect(const arrow::dataset::FileSource& source) const
{
  arrow::Schema schema{{}};
  auto fs = std::dynamic_pointer_cast<VirtualRootFileSystemBase>(source.filesystem());
  // Actually get the TTree from the ROOT file.
  auto treeFs = std::dynamic_pointer_cast<TTreeFileSystem>(fs->GetSubFilesystem(source));
  if (!treeFs.get()) {
    throw runtime_error_f("Unknown filesystem %s\n", source.filesystem()->type_name().c_str());
  }
  auto& tree = treeFs->GetTree(source);

  auto branches = tree->GetListOfBranches();
  auto n = branches->GetEntries();

  std::vector<BranchInfo> branchInfos;
  for (auto i = 0; i < n; ++i) {
    auto branch = static_cast<TBranch*>(branches->At(i));
    auto name = std::string{branch->GetName()};
    auto pos = name.find("_size");
    if (pos != std::string::npos) {
      name.erase(pos);
      branchInfos.emplace_back(BranchInfo{name, (TBranch*)nullptr, true});
    } else {
      auto lookup = std::find_if(branchInfos.begin(), branchInfos.end(), [&](BranchInfo const& bi) {
        return bi.name == name;
      });
      if (lookup == branchInfos.end()) {
        branchInfos.emplace_back(BranchInfo{name, branch, false});
      } else {
        lookup->ptr = branch;
      }
    }
  }

  std::vector<std::shared_ptr<arrow::Field>> fields;
  tree->SetCacheSize(25000000);
  for (auto& bi : branchInfos) {
    static TClass* cls;
    EDataType type;
    bi.ptr->GetExpectedType(cls, type);
    auto listSize = -1;
    if (!bi.mVLA) {
      listSize = static_cast<TLeaf*>(bi.ptr->GetListOfLeaves()->At(0))->GetLenStatic();
    }
    auto field = std::make_shared<arrow::Field>(bi.ptr->GetName(), arrowTypeFromROOT(type, listSize));
    fields.push_back(field);

    tree->AddBranchToCache(bi.ptr);
    if (strncmp(bi.ptr->GetName(), "fIndexArray", strlen("fIndexArray")) == 0) {
      std::string sizeBranchName = bi.ptr->GetName();
      sizeBranchName += "_size";
      auto* sizeBranch = (TBranch*)tree->GetBranch(sizeBranchName.c_str());
      if (sizeBranch) {
        tree->AddBranchToCache(sizeBranch);
      }
    }
  }
  tree->StopCacheLearningPhase();

  return std::make_shared<arrow::Schema>(fields);
}

/// \brief Create a FileFragment for a FileSource.
arrow::Result<std::shared_ptr<arrow::dataset::FileFragment>> TTreeFileFormat::MakeFragment(
  arrow::dataset::FileSource source, arrow::compute::Expression partition_expression,
  std::shared_ptr<arrow::Schema> physical_schema)
{
  std::shared_ptr<arrow::dataset::FileFormat> format = std::make_shared<TTreeFileFormat>(mTotCompressedSize, mTotUncompressedSize);

  auto fragment = std::make_shared<TTreeFileFragment>(std::move(source), std::move(format),
                                                      std::move(partition_expression),
                                                      std::move(physical_schema));
  return std::dynamic_pointer_cast<arrow::dataset::FileFragment>(fragment);
}

class TTreeFileWriter : public arrow::dataset::FileWriter
{
  std::vector<TBranch*> branches;
  std::vector<TBranch*> sizesBranches;
  std::vector<std::shared_ptr<arrow::Array>> valueArrays;
  std::vector<std::shared_ptr<arrow::Array>> sizeArrays;
  std::vector<std::shared_ptr<arrow::DataType>> valueTypes;

  std::vector<int64_t> valuesIdealBasketSize;
  std::vector<int64_t> sizeIdealBasketSize;

  std::vector<int64_t> typeSizes;
  std::vector<int64_t> listSizes;
  bool firstBasket = true;

  // This is to create a batsket size according to the first batch.
  void finaliseBasketSize(std::shared_ptr<arrow::RecordBatch> firstBatch)
  {
    O2_SIGNPOST_ID_FROM_POINTER(sid, root_arrow_fs, this);
    O2_SIGNPOST_START(root_arrow_fs, sid, "finaliseBasketSize", "First batch with %lli rows received and %zu columns",
                      firstBatch->num_rows(), firstBatch->columns().size());
    for (size_t i = 0; i < branches.size(); i++) {
      auto* branch = branches[i];
      auto* sizeBranch = sizesBranches[i];

      int valueSize = valueTypes[i]->byte_width();
      if (listSizes[i] == 1) {
        O2_SIGNPOST_EVENT_EMIT(root_arrow_fs, sid, "finaliseBasketSize", "Branch %s exists and uses %d bytes per entry for %lli entries.",
                               branch->GetName(), valueSize, firstBatch->num_rows());
        assert(sizeBranch == nullptr);
        branch->SetBasketSize(1024 + firstBatch->num_rows() * valueSize);
      } else if (listSizes[i] == -1) {
        O2_SIGNPOST_EVENT_EMIT(root_arrow_fs, sid, "finaliseBasketSize", "Branch %s exists and uses %d bytes per entry.",
                               branch->GetName(), valueSize);
        // This should probably lookup the
        auto column = firstBatch->GetColumnByName(schema_->field(i)->name());
        auto list = std::static_pointer_cast<arrow::ListArray>(column);
        O2_SIGNPOST_EVENT_EMIT(root_arrow_fs, sid, "finaliseBasketSize", "Branch %s needed. Associated size branch %s and there are %lli entries of size %d in that list.",
                               branch->GetName(), sizeBranch->GetName(), list->length(), valueSize);
        branch->SetBasketSize(1024 + firstBatch->num_rows() * valueSize * list->length());
        sizeBranch->SetBasketSize(1024 + firstBatch->num_rows() * 4);
      } else {
        O2_SIGNPOST_EVENT_EMIT(root_arrow_fs, sid, "finaliseBasketSize", "Branch %s needed. There are %lli entries per array of size %d in that list.",
                               branch->GetName(), listSizes[i], valueSize);
        assert(sizeBranch == nullptr);
        branch->SetBasketSize(1024 + firstBatch->num_rows() * valueSize * listSizes[i]);
      }

      auto field = firstBatch->schema()->field(i);
      if (field->name().starts_with("fIndexArray")) {
        // One int per array to keep track of the size
        int idealBasketSize = 4 * firstBatch->num_rows() + 1024 + field->type()->byte_width() * firstBatch->num_rows(); // minimal additional size needed, otherwise we get 2 baskets
        int basketSize = std::max(32000, idealBasketSize);                                                              // keep a minimum value
        sizeBranch->SetBasketSize(basketSize);
        branch->SetBasketSize(basketSize);
      }
    }
    O2_SIGNPOST_END(root_arrow_fs, sid, "finaliseBasketSize", "Done");
  }

 public:
  // Create the TTree based on the physical_schema, not the one in the batch.
  // The write method will have to reconcile the two schemas.
  TTreeFileWriter(std::shared_ptr<arrow::Schema> schema, std::shared_ptr<arrow::dataset::FileWriteOptions> options,
                  std::shared_ptr<arrow::io::OutputStream> destination,
                  arrow::fs::FileLocator destination_locator)
    : FileWriter(schema, options, destination, destination_locator)
  {
    // Batches have the same number of entries for each column.
    auto directoryStream = std::dynamic_pointer_cast<TDirectoryFileOutputStream>(destination_);
    auto treeStream = std::dynamic_pointer_cast<TTreeOutputStream>(destination_);

    if (directoryStream.get()) {
      TDirectoryFile* dir = directoryStream->GetDirectory();
      dir->cd();
      auto* tree = new TTree(destination_locator_.path.c_str(), "");
      treeStream = std::make_shared<TTreeOutputStream>(tree, "");
    } else if (treeStream.get()) {
      // We already have a tree stream, let's derive a new one
      // with the destination_locator_.path as prefix for the branches
      // This way we can multiplex multiple tables in the same tree.
      auto* tree = treeStream->GetTree();
      treeStream = std::make_shared<TTreeOutputStream>(tree, destination_locator_.path);
    } else {
      // I could simply set a prefix here to merge to an already existing tree.
      throw std::runtime_error("Unsupported backend.");
    }

    for (auto i = 0u; i < schema->fields().size(); ++i) {
      auto& field = schema->field(i);
      listSizes.push_back(1);

      int valuesIdealBasketSize = 0;
      // Construct all the needed branches.
      switch (field->type()->id()) {
        case arrow::Type::FIXED_SIZE_LIST: {
          listSizes.back() = std::static_pointer_cast<arrow::FixedSizeListType>(field->type())->list_size();
          valuesIdealBasketSize = 1024 + valueTypes.back()->byte_width() * listSizes.back();
          valueTypes.push_back(field->type()->field(0)->type());
          sizesBranches.push_back(nullptr);
          std::string leafList = fmt::format("{}[{}]{}", field->name(), listSizes.back(), rootSuffixFromArrow(valueTypes.back()->id()));
          branches.push_back(treeStream->CreateBranch(field->name().c_str(), leafList.c_str()));
        } break;
        case arrow::Type::LIST: {
          valueTypes.push_back(field->type()->field(0)->type());
          std::string leafList = fmt::format("{}[{}_size]{}", field->name(), field->name(), rootSuffixFromArrow(valueTypes.back()->id()));
          listSizes.back() = -1; // VLA, we need to calculate it on the fly;
          std::string sizeLeafList = field->name() + "_size/I";
          sizesBranches.push_back(treeStream->CreateBranch((field->name() + "_size").c_str(), sizeLeafList.c_str()));
          branches.push_back(treeStream->CreateBranch(field->name().c_str(), leafList.c_str()));
          // Notice that this could be replaced by a better guess of the
          // average size of the list elements, but this is not trivial.
        } break;
        default: {
          valueTypes.push_back(field->type());
          std::string leafList = field->name() + rootSuffixFromArrow(valueTypes.back()->id());
          sizesBranches.push_back(nullptr);
          branches.push_back(treeStream->CreateBranch(field->name().c_str(), leafList.c_str()));
        } break;
      }
    }
    // We create the branches from the schema
  }

  arrow::Status Write(const std::shared_ptr<arrow::RecordBatch>& batch) override
  {
    if (firstBasket) {
      firstBasket = false;
      finaliseBasketSize(batch);
    }

    // Support writing empty tables
    if (batch->columns().empty() || batch->num_rows() == 0) {
      return arrow::Status::OK();
    }

    // Batches have the same number of entries for each column.
    auto directoryStream = std::dynamic_pointer_cast<TDirectoryFileOutputStream>(destination_);
    TTree* tree = nullptr;
    if (directoryStream.get()) {
      TDirectoryFile* dir = directoryStream->GetDirectory();
      tree = (TTree*)dir->Get(destination_locator_.path.c_str());
    }
    auto treeStream = std::dynamic_pointer_cast<TTreeOutputStream>(destination_);

    if (!tree) {
      // I could simply set a prefix here to merge to an already existing tree.
      throw std::runtime_error("Unsupported backend.");
    }

    for (auto i = 0u; i < batch->columns().size(); ++i) {
      auto column = batch->column(i);
      auto& field = batch->schema()->field(i);

      valueArrays.push_back(nullptr);

      switch (field->type()->id()) {
        case arrow::Type::FIXED_SIZE_LIST: {
          auto list = std::static_pointer_cast<arrow::FixedSizeListArray>(column);
          valueArrays.back() = list->values();
        } break;
        case arrow::Type::LIST: {
          auto list = std::static_pointer_cast<arrow::ListArray>(column);
          valueArrays.back() = list;
        } break;
        case arrow::Type::BOOL: {
          // In case of arrays of booleans, we need to go back to their
          // char based representation for ROOT to save them.
          auto boolArray = std::static_pointer_cast<arrow::BooleanArray>(column);

          int64_t length = boolArray->length();
          arrow::UInt8Builder builder;
          auto ok = builder.Reserve(length);

          for (int64_t i = 0; i < length; ++i) {
            if (boolArray->IsValid(i)) {
              // Expand each boolean value (true/false) to uint8 (1/0)
              uint8_t value = boolArray->Value(i) ? 1 : 0;
              auto ok = builder.Append(value);
            } else {
              // Append null for invalid entries
              auto ok = builder.AppendNull();
            }
          }
          valueArrays.back() = *builder.Finish();
        } break;
        default:
          valueArrays.back() = column;
      }
    }

    int64_t pos = 0;
    while (pos < batch->num_rows()) {
      for (size_t bi = 0; bi < branches.size(); ++bi) {
        auto* branch = branches[bi];
        auto* sizeBranch = sizesBranches[bi];
        auto array = batch->column(bi);
        auto& field = batch->schema()->field(bi);
        auto& listSize = listSizes[bi];
        auto valueType = valueTypes[bi];
        auto valueArray = valueArrays[bi];

        switch (field->type()->id()) {
          case arrow::Type::LIST: {
            auto list = std::static_pointer_cast<arrow::ListArray>(array);
            listSize = list->value_length(pos);
            uint8_t const* buffer = std::static_pointer_cast<arrow::PrimitiveArray>(valueArray)->values()->data() + array->offset() + list->value_offset(pos) * valueType->byte_width();
            branch->SetAddress((void*)buffer);
            sizeBranch->SetAddress(&listSize);
          };
            break;
          case arrow::Type::FIXED_SIZE_LIST:
          default: {
            uint8_t const* buffer = std::static_pointer_cast<arrow::PrimitiveArray>(valueArray)->values()->data() + array->offset() + pos * listSize * valueType->byte_width();
            branch->SetAddress((void*)buffer);
          };
        }
      }
      tree->Fill();
      ++pos;
    }
    return arrow::Status::OK();
  }

  arrow::Future<> FinishInternal() override
  {
    auto treeStream = std::dynamic_pointer_cast<TTreeOutputStream>(destination_);
    auto* tree = treeStream->GetTree();
    tree->Write("", TObject::kOverwrite);
    tree->SetDirectory(nullptr);

    return {};
  };
};
arrow::Result<std::shared_ptr<arrow::dataset::FileWriter>> TTreeFileFormat::MakeWriter(std::shared_ptr<arrow::io::OutputStream> destination, std::shared_ptr<arrow::Schema> schema, std::shared_ptr<arrow::dataset::FileWriteOptions> options, arrow::fs::FileLocator destination_locator) const
{
  auto writer = std::make_shared<TTreeFileWriter>(schema, options, destination, destination_locator);
  return std::dynamic_pointer_cast<arrow::dataset::FileWriter>(writer);
}

std::shared_ptr<arrow::dataset::FileWriteOptions> TTreeFileFormat::DefaultWriteOptions()
{
  std::shared_ptr<TTreeFileWriteOptions> options(
    new TTreeFileWriteOptions(shared_from_this()));
  return options;
}

TTreeFileSystem::~TTreeFileSystem() = default;

DEFINE_DPL_PLUGINS_BEGIN
DEFINE_DPL_PLUGIN_INSTANCE(TTreeObjectReadingImplementation, RootObjectReadingImplementation);
DEFINE_DPL_PLUGINS_END
} // namespace o2::framework
