#include <mpi.h>
#include <cstdio>
#include <cstdarg>
#include <cstring>
#include "mpi_logging.h"

struct mpi_log_file_s {
	MPI_File handle;
	bool shared;
};

static char FilenamePrefix[128];
static int FilenamePrefixLen;
static bool HasFilenamePrefix;

static MPI_File DefaultLogFile;

static int CurrentRank;
static char CurrentRankStr[16];
static int CurrentRankStrLen;

static void generate_full_filename(const char* filename_in, char* full_filename_out, int full_filename_size, bool with_rank)
{
	char* out_end = full_filename_out + full_filename_size - 1;
	char* out = full_filename_out;
	*out = '\0';
	if (HasFilenamePrefix) {
		strncpy(out, FilenamePrefix, out_end-out);
		out = full_filename_out + strlen(full_filename_out);
		if(out>=out_end) return;
		*out = '.';
		out++;
		if(out>=out_end) return;
	}
	if (with_rank) {
		strncpy(out, CurrentRankStr, out_end-out);
		out = full_filename_out + strlen(full_filename_out);
		if(out>=out_end) return;
		*out = '.';
		out++;
		if(out>=out_end) return;
	}
	strncpy(out, filename_in, out_end-out);
	*out_end = '\0';
}

void init_mpi_logging(const char* default_log_filename, const char* filename_prefix)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &CurrentRank);
	snprintf(CurrentRankStr, 16, "%d", CurrentRank);
	CurrentRankStr[15] = 0;
	CurrentRankStrLen = strlen(CurrentRankStr);

	strncpy(FilenamePrefix, filename_prefix, 128);
	FilenamePrefix[127] = 0;
	FilenamePrefixLen = strlen(FilenamePrefix);
	if(FilenamePrefixLen==0) HasFilenamePrefix=false;
	else HasFilenamePrefix=true;

	MPI_File_open(MPI_COMM_WORLD, default_log_filename, MPI_MODE_APPEND | MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &DefaultLogFile);
}

static mpi_log_file mpi_open_log_file_internal(const char* filename, bool shared)
{
	char full_filename[256];
	MPI_File handle;
	int ret;
	generate_full_filename(filename, full_filename, 256, !shared);
	ret = MPI_File_open(MPI_COMM_WORLD, full_filename, MPI_MODE_APPEND | MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &handle);
	if(ret != MPI_SUCCESS) return nullptr;

	mpi_log_file_s* f = new mpi_log_file_s;
	f->handle = handle;
	f->shared = shared;
	return f;
}

mpi_log_file mpi_open_log_file(const char* filename)
{
	return mpi_open_log_file_internal(filename, false);
}

mpi_log_file mpi_open_log_file_shared(const char* filename)
{
	return mpi_open_log_file_internal(filename, true);
}

void mpi_close_log_file(mpi_log_file f)
{
	MPI_File_close(f->handle);
	delete f;
}

void mpi_log(mpi_log_file file, const char *fmt, ...)
{
	bool defaultLog = !file;
	char str[512];
	int start_index = 0;
	int count;
	va_list args;

	if (defaultLog) {
		snprintf(str, 511, "[Rank %d] ", CurrentRank);
		start_index = strlen(str);
	}

	va_start(args, fmt);
	count = vsnprintf(str+start_index, 511-start_index, fmt, args);
	start_index += count;
	if(start_index>510) start_index=510;
	va_end(args);
	str[start_index] = '\n';
	start_index++;
	str[start_index] = 0;

	if (defaultLog) {
		MPI_File_write_shared(DefaultLogFile, str, start_index, MPI_CHAR, MPI_STATUS_IGNORE);
	} else if (file->shared) {
		MPI_File_write_shared(file->handle, str, start_index, MPI_CHAR, MPI_STATUS_IGNORE);
	} else {
		MPI_File_write(file->handle, str, start_index, MPI_CHAR, MPI_STATUS_IGNORE);
	}
}
