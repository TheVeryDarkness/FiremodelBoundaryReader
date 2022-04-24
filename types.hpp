#pragma once
#include <exception>
#include <vector>

using std::exception;
using std::vector;
class file_error : public std::exception {
public:
  const char *what() const noexcept { return "File is not correct. Check it."; }
};
class stream_error : public std::exception {
public:
  const char *what() const noexcept {
    return "File stream terminated unexpectedly. This is usually because of an "
           "unexpected EOF.";
  }
};

#define CHECK_STREAM(STREAM)                                                   \
  if (!(STREAM))                                                               \
    throw stream_error();

#define CHECK_FORMAT(EXPR)                                                     \
  if (!(EXPR))                                                                 \
    throw file_error();

using i64 = std::int64_t;
using u64 = std::uint64_t;
using i32 = std::int32_t;
using u32 = std::uint32_t;
using u16 = std::uint16_t;