#pragma once

#include <sys/syscall.h>
#include <sys/types.h>
#include <unistd.h>

#include <cassert>
#include <iostream>
#define gettid() syscall(SYS_gettid)

#define LOG_FATAL(log_fmt, log_arg...)                                   \
  do {                                                                   \
    std::printf("Fatal [%s:%d][%s][thread %ld] " log_fmt "\n", __FILE__, \
                __LINE__, __FUNCTION__, gettid(), ##log_arg);            \
    assert(0);                                                           \
  } while (0)

#define LOG_ERROR(log_fmt, log_arg...)                                   \
  do {                                                                   \
    std::printf("Error [%s:%d][%s][thread %ld] " log_fmt "\n", __FILE__, \
                __LINE__, __FUNCTION__, gettid(), ##log_arg);            \
  } while (0)

#define LOG_WARN(log_fmt, log_arg...)                                      \
  do {                                                                     \
    std::printf("Warning [%s:%d][%s][thread %ld] " log_fmt "\n", __FILE__, \
                __LINE__, __FUNCTION__, gettid(), ##log_arg);              \
  } while (0)

#define LOG_INFO(log_fmt, log_arg...)                                   \
  do {                                                                  \
    std::printf("Info [%s:%d][%s][thread %ld] " log_fmt "\n", __FILE__, \
                __LINE__, __FUNCTION__, gettid(), ##log_arg);           \
  } while (0)

// #define CASE_DEBUG