#pragma once

typedef struct mpi_log_file_s *mpi_log_file;

void init_mpi_logging(const char* default_log_filename, const char* filename_prefix); // Pass "" for no prefix. Prefix is not applied to the default log file. Prefix is limited to 127 characters.
mpi_log_file mpi_open_log_file(const char* filename); // File name is <filename_prefix>.<Current MPI Rank>.<filename>. The total filename is limited to 255 characters.
mpi_log_file mpi_open_log_file_shared(const char* filename); // File name is <filename_prefix>.<filename>. The total filename is limited to 255 characters.
void mpi_close_log_file(mpi_log_file file);

// If file is NULL, msg is logged to the default log file. A newline character is automatically appended.
// Msgs logged to the default log file have "[Rank %d] " prepended, with %d replaced with the current MPI rank.
// The logged msg is truncated to 511 characters including the newline character and the length of the "[Rank %d]" string, if present.
void mpi_log(mpi_log_file file, const char* fmt, ...);
