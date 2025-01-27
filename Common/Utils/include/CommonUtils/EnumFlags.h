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
#ifndef O2_FRAMEWORK_FLAGS_H_
#define O2_FRAMEWORK_FLAGS_H_

#include <algorithm>
#include <array>
#include <concepts>
#include <exception>
#include <ostream>
#include <source_location>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <string>
#include <sstream>
#include <limits>
#include <bitset>
#include <initializer_list>
#include <cstdint>
#include <cstddef>
#include <cctype>
#include <utility>
#include <optional>
#include <iostream>
#include <iomanip>

#include "CommonUtils/StringUtils.h"

namespace o2::utils
{

namespace details::enum_flags
{

// Require that an enum with an underlying unsigned type.
template <typename E>
concept EnumFlagHelper = requires {
  requires std::is_enum_v<E>;
  requires std::is_unsigned_v<std::underlying_type_t<E>>;
  requires std::same_as<E, std::decay_t<E>>;
};

// Static constexpr only helper struct to implement modicum of enum reflection
// functions and also check via concepts expected properties of the enum.
// This is very much inspired by much more extensive libraries like magic_enum.
// Inspiration by its c++20 version (https://github.com/fix8mt/conjure_enum).
template <EnumFlagHelper E>
struct FlagsHelper final {
  using U = std::underlying_type_t<E>;

  static constexpr bool isScoped() noexcept
  {
    return std::is_enum_v<E> && !std::is_convertible_v<E, std::underlying_type_t<E>>;
  }

  // Return line at given position.
  template <E e>
  static consteval const char* tpeek() noexcept
  {
    return std::source_location::current().function_name();
  }
  // string_view value of function above
  template <E e>
  static constexpr std::string_view tpeek_v{tpeek<e>()};

  // Compiler Specifics
  static constexpr auto CSpecifics{std::to_array<
    std::tuple<std::string_view, char, std::string_view, char>>({
#if defined __clang__
    {"e = ", ']', "(anonymous namespace)", '('},
    {"T = ", ']', "(anonymous namespace)", '('},
#else // assuming __GNUC__
    {"e = ", ';', "<unnamed>", '<'},
    {"T = ", ']', "{anonymous}", '{'},
#endif
  })};
  enum class SVal : uint8_t { Start,
                              End,
                              AnonStr,
                              AnonStart };
  enum class SType : uint8_t { Enum_t,
                               Type_t,
                               eT0,
                               eT1,
                               eT2,
                               eT3 };
  // Extract a compiler specification.
  template <SVal v, SType t>
  static constexpr auto getSpec() noexcept
  {
    return std::get<static_cast<size_t>(v)>(CSpecifics[static_cast<size_t>(t)]);
  }

  // Range that is scanned by the compiler
  static constexpr size_t MinScan{0};
  static constexpr size_t MarginScan{1};                                // Scan one past to check for overpopulation
  static constexpr size_t MaxUnderScan{std::numeric_limits<U>::digits}; // Maximum digits the underlying type has
  static constexpr size_t MaxScan{MaxUnderScan + MarginScan};

  // Checks if a given 'localation' contains an enum.
  template <E e>
  static constexpr bool isValid() noexcept
  {
    constexpr auto tp{tpeek_v<e>.rfind(getSpec<SVal::Start, SType::Enum_t>())};
    if constexpr (tp == std::string_view::npos) {
      return false;
    }
#if defined __clang__
    else if constexpr (tpeek_v<e>[tp + getSpec<SVal::Start, SType::Enum_t>().size()] == '(') {
      if constexpr (tpeek_v<e>[tp + getSpec<SVal::Start, SType::Enum_t>().size() + 1] == '(') {
        return false;
      }
      if constexpr (tpeek_v<e>.find(getSpec<SVal::AnonStr, SType::Enum_t>(), tp + getSpec<SVal::Start, SType::Enum_t>().size()) != std::string_view::npos) {
        return true;
      }
    } else if constexpr (tpeek_v<e>.find_first_of(getSpec<SVal::End, SType::Enum_t>(), tp + getSpec<SVal::Start, SType::Enum_t>().size()) != std::string_view::npos) {
      // check if this is an anonymous enum
      return true;
    }
    return false;
#else
    else if constexpr (tpeek_v<e>[tp + getSpec<SVal::Start, SType::Enum_t>().size()] != '(' && tpeek_v<e>.find_first_of(getSpec<SVal::End, SType::Enum_t>(), tp + getSpec<SVal::Start, SType::Enum_t>().size()) != std::string_view::npos) {
      return true;
    } else {
      return false;
    }
#endif
  }

  // Extract which values are present in the enum by checking all values in
  // the min-max-range above.
  template <size_t... I>
  static constexpr auto getValues(std::index_sequence<I...> /*unused*/) noexcept
  {
    constexpr std::array<bool, sizeof...(I)> valid{isValid<static_cast<E>(MinScan + I)>()...};
    constexpr auto count{std::count_if(valid.cbegin(), valid.cend(), [](bool v) noexcept { return v; })};
    static_assert(count > 0, "Requiring non-empty enum!");
    static_assert(count <= MaxUnderScan, "Underlying type of enum has less digits than given expected!");
    std::array<E, count> values{};
    for (size_t idx{}, n{}; n < count; ++idx) {
      if (valid[idx]) {
        values[n++] = static_cast<E>(MinScan + idx);
      }
    }
    return values;
  }
  static constexpr auto Values{getValues(std::make_index_sequence<MaxScan - MinScan - MarginScan>())};              // Enum Values
  static constexpr auto count() noexcept { return Values.size(); }                                                  // Number of enum members
  static constexpr auto Min_v{Values.front()};                                                                      // Enum first entry
  static constexpr auto Max_v{Values.back()};                                                                       // Enum last entry
  static constexpr auto Min_u_v{static_cast<size_t>(Min_v)};                                                        // Enum first entry as size_t
  static constexpr auto Max_u_v{static_cast<size_t>(Max_v)};                                                        // Enum last entry as size_t
  static constexpr bool isContinuous() noexcept { return (Max_u_v - Min_u_v + 1) == count(); }                      // Is the enum continuous
  static constexpr uint64_t MaxRep{(Max_u_v >= 64) ? std::numeric_limits<uint64_t>::max() : (1ULL << Max_u_v) - 1}; // largest representable value

  template <E e>
  static constexpr std::string_view getName()
  {
    constexpr auto tp{tpeek_v<e>.rfind(getSpec<SVal::Start, SType::Enum_t>())};
    if constexpr (tp == std::string_view::npos) {
      return {};
    }
    if constexpr (tpeek_v<e>[tp + getSpec<SVal::Start, SType::Enum_t>().size()] == getSpec<SVal::AnonStart, SType::Enum_t>()) {
#if defined __clang__
      if constexpr (tpeek_v<e>[tp + getSpec<SVal::Start, SType::enum_t>().size() + 1] == getSpec<SVal::AnonStart, SType::Enum_t>()) {
        return {};
      }
#endif
      if (constexpr auto lstr{tpeek_v<e>.substr(tp + getSpec<SVal::Start, SType::Enum_t>().size())}; lstr.find(getSpec<SVal::AnonStr, SType::Enum_t>()) != std::string_view::npos) { // is anon
        if constexpr (constexpr auto lc{lstr.find_first_of(getSpec<SVal::End, SType::Enum_t>())}; lc != std::string_view::npos) {
          return lstr.substr(getSpec<SVal::AnonStr, SType::Enum_t>().size() + 2, lc - (getSpec<SVal::AnonStr, SType::Enum_t>().size() + 2));
        }
      }
    }
    constexpr std::string_view result{tpeek_v<e>.substr(tp + getSpec<SVal::Start, SType::Enum_t>().size())};
    if constexpr (constexpr auto lc{result.find_first_of(getSpec<SVal::End, SType::Enum_t>())}; lc != std::string_view::npos) {
      return result.substr(0, lc);
    } else {
      return {};
    }
  }

  static constexpr std::string_view removeScope(std::string_view s)
  {
    if (const auto lc{s.find_last_of(':')}; lc != std::string_view::npos) {
      return s.substr(lc + 1);
    }
    return s;
  }

  static constexpr std::string_view findScope(std::string_view s)
  {
    const auto pos1 = s.rfind("::");
    if (pos1 == std::string_view::npos) {
      return s;
    }
    const auto pos2 = s.rfind("::", pos1 - 1);
    if (pos2 == std::string_view::npos) {
      return s.substr(0, pos1);
    }
    return s.substr(pos2 + 2, pos1 - pos2 - 2);
  }

  template <E e>
  static constexpr auto getNameValue{getName<e>()};

  template <bool with_scope, std::size_t... I>
  static constexpr auto getNames(std::index_sequence<I...> /*unused*/)
  {
    if constexpr (with_scope) {
      return std::array<std::string_view, sizeof...(I)>{getNameValue<Values[I]>...};
    } else {
      return std::array<std::string_view, sizeof...(I)>{removeScope(getNameValue<Values[I]>)...};
    }
  }

  static constexpr auto Names{getNames<false>(std::make_index_sequence<count()>())};      // Enum names without scope
  static constexpr auto NamesScoped{getNames<true>(std::make_index_sequence<count()>())}; // Enum names with scope
  static constexpr auto Scope{findScope(NamesScoped.front())};                            // Enum scope

  static constexpr auto getLongestName() noexcept
  {
    size_t max{0};
    for (size_t i{0}; i < count(); ++i) {
      max = std::max(max, Names[i].size());
    }
    return max;
  }

  static constexpr auto NamesLongest{getLongestName()}; // Size of longest name

  template <E e>
  static constexpr std::string_view toString() noexcept
  {
    return getNameValue<e>();
  }

  static constexpr std::optional<E> fromString(std::string_view str) noexcept
  {
    for (std::size_t i{0}; i < count(); ++i) {
      if (Names[i] == str || NamesScoped[i] == str) {
        return Values[i];
      }
    }
    return std::nullopt;
  }

  // Convert char to lower.
  static constexpr unsigned char toLower(const unsigned char c) noexcept
  {
    return (c >= 'A' && c <= 'Z') ? (c - 'A' + 'a') : c;
  }

  // Are these chars equal (case-insensitive).
  static constexpr bool isIEqual(const unsigned char a, const unsigned char b) noexcept
  {
    return toLower(a) == toLower(b);
  }

  // Case-insensitive comparision for string_view.
  static constexpr bool isIEqual(std::string_view s1, std::string_view s2) noexcept
  {
    if (s1.size() != s2.size()) {
      return false;
    }
    for (size_t i{0}; i < s1.size(); ++i) {
      if (!isIEqual(s1[i], s2[i])) {
        return false;
      }
    }
    return true;
  }

  static constexpr std::string_view None{"none"};
  static constexpr bool hasNone() noexcept
  {
    // check that enum does not contain memeber named 'none'
    for (size_t i{0}; i < count(); ++i) {
      if (isIEqual(Names[i], None)) {
        return true;
      }
    }
    return false;
  }

  static constexpr std::string_view All{"all"};
  static constexpr bool hasAll() noexcept
  {
    // check that enum does not contain memeber named 'all'
    for (size_t i{0}; i < count(); ++i) {
      if (isIEqual(Names[i], All)) {
        return true;
      }
    }
    return false;
  }
};

} // namespace details::enum_flags

// Require an enum to fullfil what one would except from a bitset.
template <typename E>
concept EnumFlag = requires {
  // range checks
  requires details::enum_flags::FlagsHelper<E>::Min_u_v == 0;                                           // the first bit should be at position 0
  requires details::enum_flags::FlagsHelper<E>::Max_u_v < details::enum_flags::FlagsHelper<E>::count(); //  the maximum is less than the total
  requires details::enum_flags::FlagsHelper<E>::isContinuous();                                         // do not allow missing bits

  // type checks
  requires !details::enum_flags::FlagsHelper<E>::hasNone(); // added automatically
  requires !details::enum_flags::FlagsHelper<E>::hasAll();  // added automatically
};

/**
 * \brief Classs to aggregate and manage enum-based on-off flags.
 *
 * This class manages flags as bits in the underlying type of an enum, allowing
 * manipulation via enum member names. It supports operations akin to std::bitset
 * but is fully constexpr and is ideal for aggregating multiple on-off booleans,
 * e.g., enabling/disabling algorithm features.
 *
 * Example:
 * enum class AlgoOptions {
 *     Feature1,
 *     Feature2,
 *     Feature3,
 * };
 * ...
 * EnumFlags<AlgoOptions> opts;
 * opts.set("Feature1 | Feature3"); // Set Feature1 and Feature3.
 * if (opts[AlgoOptions::Feature1]) { // Do some work. } // Check if Feature1 is set.
 *
 * Additional examples of how to use this class are in testEnumFlags.cxx.
 */
template <EnumFlag E>
class EnumFlags
{
  using H = details::enum_flags::FlagsHelper<E>;
  using U = std::underlying_type_t<E>;
  U mBits{0};

  // Converts enum to its underlying type.
  constexpr auto to_underlying(E e) const noexcept
  {
    return static_cast<U>(e);
  }

  // Returns the bit representation of a flag.
  constexpr auto to_bit(E e) const noexcept
  {
    return U(1) << to_underlying(e);
  }

 public:
  // Default constructor.
  constexpr explicit EnumFlags() = default;
  // Constructor to initialize with a single flag.
  constexpr explicit EnumFlags(E e) : mBits(to_bit(e)) {}
  // Copy constructor.
  constexpr EnumFlags(const EnumFlags&) = default;
  // Move constructor.
  constexpr EnumFlags(EnumFlags&&) = default;
  // Constructor to initialize with the underlyiny type.
  constexpr explicit EnumFlags(U u) : mBits(u) {}
  // Initialize with a list of flags.
  constexpr EnumFlags(std::initializer_list<E> flags) noexcept
  {
    std::for_each(flags.begin(), flags.end(), [this](const E f) noexcept { mBits |= to_bit(f); });
  }
  // Destructor.
  constexpr ~EnumFlags() = default;

  static constexpr U None{0};        // Represents no flags set.
  static constexpr U All{H::MaxRep}; // Represents all flags set.

  // Return list of all enum values
  static constexpr auto getValues() noexcept
  {
    return H::Values;
  }

  // Return list of all enum Names
  static constexpr auto getNames() noexcept
  {
    return H::Names;
  }

  // Sets flags from a string representation.
  // This can be either from a number representation (binary or digits) or
  // a concatenation of the enums members name e.g., 'Enum1|Enum2|...'
  void set(const std::string& s, int base = 2)
  {
    // on throw restore previous state and rethrow
    const U prev = mBits;
    reset();
    try {
      setImpl(s, base);
    } catch (const std::exception& e) {
      mBits = prev;
      throw;
    }
  }
  // Returns the raw bitset value.
  constexpr auto value() const noexcept
  {
    return mBits;
  }

  // Resets all flags.
  constexpr void reset() noexcept
  {
    mBits = U(0);
  }

  // Resets a specific flag.
  template <typename T>
    requires std::is_same_v<T, E>
  constexpr void reset(T t)
  {
    mBits &= ~to_bit(t);
  }

  // Tests if a specific flag is set.
  template <typename T>
    requires std::is_same_v<T, E>
  [[nodiscard]] constexpr bool test(T t) const noexcept
  {
    return (mBits & to_bit(t)) != None;
  }

  // Sets a specific flag.
  template <typename T>
    requires std::is_same_v<T, E>
  constexpr void set(T t) noexcept
  {
    mBits |= to_bit(t);
  }

  // Toggles a specific flag.
  template <typename T>
    requires std::is_same_v<T, E>
  constexpr void toggle(T t) noexcept
  {
    mBits ^= to_bit(t);
  }

  // Checks if any flag is set.
  [[nodiscard]] constexpr bool any() const noexcept
  {
    return mBits != None;
  }

  // Returns the bitset as a binary string.
  [[nodiscard]] std::string string() const
  {
    std::ostringstream oss;
    oss << std::bitset<H::count()>(mBits);
    return oss.str();
  }

  // Returns the bitset as a pretty multiline binary string.
  [[nodiscard]] std::string pstring(bool withNewline = false) const
  {
    std::ostringstream oss;
    if (withNewline) {
      oss << '\n';
    }
    oss << "0b";
    const std::bitset<H::count()> bits(mBits);
    oss << bits;
    if constexpr (H::isScoped()) {
      oss << " " << H::Scope;
    }
    oss << '\n';
    for (size_t i = 0; i < H::count(); ++i) {
      oss << "  ";
      for (size_t j = 0; j < H::count() - i - 1; ++j) {
        oss << "┃";
      }
      oss << "┗";
      for (size_t a{2 + i}; --a != 0U;) {
        oss << "━";
      }
      oss << " " << std::setw(H::NamesLongest) << std::left
          << H::Names[i] << " " << (bits[i] ? "[Active]" : "[Inactive]");
      if (i != H::count() - 1) {
        oss << "\n";
      }
    }
    return oss.str();
  }

  // Checks if any flag is set (Boolean context).
  constexpr explicit operator bool() const noexcept
  {
    return any();
  }

  // Check if given flag is set.
  template <typename T>
    requires std::is_same_v<T, E>
  constexpr bool operator[](const T t) noexcept
  {
    return test(t);
  }

  // Checks if two flag sets are equal.
  constexpr bool operator==(const EnumFlags& o) const noexcept
  {
    return mBits == o.mBits;
  }

  // Checks if two flag sets are not equal.
  constexpr bool operator!=(const EnumFlags& o) const noexcept
  {
    return mBits != o.mBits;
  }

  // Copy assignment operator
  constexpr EnumFlags& operator=(const EnumFlags& o) = default;

  // Move assignment operator
  constexpr EnumFlags& operator=(EnumFlags&& o) = default;

  // Performs a bitwise OR with a flag.
  template <typename T>
    requires std::is_same_v<T, E>
  constexpr EnumFlags& operator|=(T t) noexcept
  {
    mBits |= to_bit(t);
    return *this;
  }

  // Performs a bitwise AND with a flag.
  template <typename T>
    requires std::is_same_v<T, E>
  constexpr EnumFlags& operator&=(T t) noexcept
  {
    mBits &= to_bit(t);
    return *this;
  }

  // Returns a flag set with a bitwise AND.
  template <typename T>
    requires std::is_same_v<T, E>
  constexpr EnumFlags operator&(T t) const noexcept
  {
    return EnumFlags(mBits & to_bit(t));
  }

  // Returns a flag set with all bits inverted.
  constexpr EnumFlags operator~() const noexcept
  {
    return EnumFlags(~mBits);
  }

  // Performs a bitwise OR with another flag set.
  constexpr EnumFlags operator|(const EnumFlags& o) const noexcept
  {
    return EnumFlags(mBits | o.mBits);
  }

  // Performs a bitwise OR assignment.
  constexpr EnumFlags& operator|=(const EnumFlags& o) noexcept
  {
    mBits |= o.mBits;
    return *this;
  }

  // Performs a bitwise XOR with another flag set.
  constexpr EnumFlags operator^(const EnumFlags& o) const noexcept
  {
    return Flags(mBits ^ o.mBits);
  }

  // Performs a bitwise XOR assignment.
  constexpr EnumFlags& operator^=(const EnumFlags& o) noexcept
  {
    mBits ^= o.mBits;
    return *this;
  }

  // Checks if all specified flags are set.
  template <typename... Ts>
  constexpr bool all_of(Ts... flags) const noexcept
  {
    return ((test(flags) && ...));
  }

  // Checks if none of the specified flags are set.
  template <typename... Ts>
  constexpr bool none_of(Ts... flags) const noexcept
  {
    return (!(test(flags) || ...));
  }

  // Serializes the flag set to a string.
  [[nodiscard]] std::string serialize() const
  {
    return std::to_string(mBits);
  }

  // Deserializes a string into the flag set.
  void deserialize(const std::string& data)
  {
    uint64_t v = std::stoul(data);
    if (v > H::MaxRep) {
      throw std::out_of_range("Values exceeds enum range.");
    }
    mBits = static_cast<U>(v);
  }

  // Counts the number of set bits (active flags).
  [[nodiscard]] constexpr size_t count() const noexcept
  {
    size_t c{0};
    for (size_t i{H::Min_u_v}; i < H::Max_u_v; ++i) {
      if ((mBits & (U(1) << i)) != U(0)) {
        ++c;
      }
    }
    return c;
  }

  // Returns the union of two flag sets.
  constexpr EnumFlags union_with(const EnumFlags& o) const noexcept
  {
    return EnumFlags(mBits | o.mBits);
  }

  // Returns the intersection of two flag sets.
  constexpr EnumFlags intersection_with(const EnumFlags& o) const noexcept
  {
    return EnumFlags(mBits & o.mBits);
  }

  // Checks if all flags in another Flags object are present in the current object.
  constexpr bool contains(const EnumFlags& other) const noexcept
  {
    return (mBits & other.mBits) == other.mBits;
  }

 private:
  // Set implemnetation, bits was zeroed before.
  void setImpl(const std::string& s, int base = 2)
  {
    if (std::all_of(s.begin(), s.end(), [](unsigned char c) { return std::isdigit(c); })) {
      if (base == 2) { // check of only 0 and 1 in string
        if (!std::all_of(s.begin(), s.end(), [](char c) { return c == '0' || c == '1'; })) {
          throw std::invalid_argument("Invalid binary string.");
        }
      }
      uint64_t v = std::stoul(s, nullptr, base);
      if (v > H::MaxRep) {
        throw std::out_of_range("Values exceeds enum range.");
      }
      mBits = static_cast<U>(v);
    } else if (std::all_of(s.begin(), s.end(), [](unsigned char c) { return std::isalnum(c) != 0 || c == '|' || c == ' ' || c == ':'; })) {
      std::string cs{s};
      std::transform(cs.begin(), cs.end(), cs.begin(), [](unsigned char c) { return std::tolower(c); });
      if (cs == H::All) {
        mBits = All;
      } else if (cs == H::None) {
        mBits = None;
      } else {
        for (const auto& tok : Str::tokenize(s, '|')) {
          if (auto e = H::fromString(tok)) {
            mBits |= to_bit(*e);
          } else {
            throw std::invalid_argument(tok + " is not a valid enum value!");
          }
        }
      }
    } else {
      throw std::invalid_argument("Cannot parse string!");
    }
  }
};

template <EnumFlag E>
std::ostream& operator<<(std::ostream& os, const EnumFlags<E>& f)
{
  os << f.pstring(true);
  return os;
}

} // namespace o2::utils

#endif
