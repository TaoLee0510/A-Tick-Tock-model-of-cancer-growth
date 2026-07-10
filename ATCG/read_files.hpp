#ifndef read_files_hpp
#define read_files_hpp

#include <string>

#include "cell_store.hpp"
#include "cell_trace.hpp"
#include "recovery_parameters.hpp"

void read_file(CellStore &cells, CellTraceStore &cell_trace, RecoveryParameters &parameters, std::string Cell_arry_file, std::string Cell_trace_arry_file, std::string Parameters);

#endif
