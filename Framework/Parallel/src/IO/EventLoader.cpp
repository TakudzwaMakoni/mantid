#include <H5Cpp.h>

#include "MantidKernel/make_unique.h"
#include "MantidParallel/IO/EventLoader.h"
#include "MantidParallel/IO/NXEventDataLoader.h"

namespace Mantid {
namespace Parallel {
namespace IO {

EventLoader::EventLoader(const std::string &filename,
                         const std::string &groupName,
                         const std::vector<std::string> &bankNames,
                         const std::vector<int32_t> &bankOffsets)
    : m_file(Kernel::make_unique<H5::H5File>(filename, H5F_ACC_RDONLY)),
      m_groupName(groupName), m_bankNames(bankNames),
      m_bankOffsets(bankOffsets) {}

EventLoader::~EventLoader() = default;

void EventLoader::load() {
  load(readDataType("event_time_zero"), readDataType("event_time_offset"));
}

void EventLoader::load(const H5::DataType &timeZeroType,
                       const H5::DataType &timeOffsetType) {
  // Translate from DataType to actual type. Done in two steps to avoid
  // combinatoric explosion. Step 1: event_time_zero.
  if (timeZeroType == H5::PredType::NATIVE_INT64)
    return load<int64_t>(timeOffsetType);
  if (timeZeroType == H5::PredType::NATIVE_DOUBLE)
    return load<double>(timeOffsetType);
  throw std::runtime_error("Unsupported H5::DataType of event_time_zero");
}

template <class TimeZeroType>
void EventLoader::load(const H5::DataType &timeOffsetType) {
  // Translate from DataType to actual type. Done in two steps to avoid
  // combinatoric explosion. Step 2: event_time_offset.
  if (timeOffsetType == H5::PredType::NATIVE_INT32)
    return load<TimeZeroType, int32_t>();
  if (timeOffsetType == H5::PredType::NATIVE_FLOAT)
    return load<TimeZeroType, float>();
  throw std::runtime_error("Unsupported H5::DataType of event_time_offset");
}

template <class TimeZeroType, class TimeOffsetType> void EventLoader::load() {
  // determine optimal chunk size
  // make chunks
  // main load loop
}

H5::DataType EventLoader::readDataType(const std::string &name) const {
  return m_file->openDataSet(m_groupName + "/" + m_bankNames.front() + "/" +
                             name).getDataType();
}

template MANTID_PARALLEL_DLL void
EventLoader::load<int64_t>(const H5::DataType &timeOffsetType);
template MANTID_PARALLEL_DLL void
EventLoader::load<double>(const H5::DataType &timeOffsetType);

template MANTID_PARALLEL_DLL void EventLoader::load<int64_t, int32_t>();
template MANTID_PARALLEL_DLL void EventLoader::load<int64_t, float>();
template MANTID_PARALLEL_DLL void EventLoader::load<double, int32_t>();
template MANTID_PARALLEL_DLL void EventLoader::load<double, float>();

} // namespace IO
} // namespace Parallel
} // namespace Mantid
