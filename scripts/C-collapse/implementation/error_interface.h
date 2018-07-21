#ifndef ERROR_H
#define ERROR_H

#include <cstring>
#include <string>
#include <sstream>


namespace error {

void fatal_error_from_c(const char* message, int message_len);
void warning_from_c(const char* message, int message_len);

} // namespace error
#endif // ERROR_H
