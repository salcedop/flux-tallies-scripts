#ifndef HDF5_INTERFACE_H
#define HDF5_INTERFACE_H


#include "hdf5.h"
#include "hdf5_hl.h"

#include <array>
#include <string>
#include <sstream>
#include <complex.h>

namespace hdf5 {

hid_t file_open(const std::string& filename, char mode);
void file_close(hid_t file_id);
void get_name(hid_t obj_id, char* name);
int get_num_datasets(hid_t group_id);
int get_num_groups(hid_t group_id);
void get_datasets(hid_t group_id, char* name[]);
void get_groups(hid_t group_id, char* name[]);
void get_shape(hid_t obj_id, hsize_t* dims);
void get_shape_attr(hid_t obj_id, const char* name, hsize_t* dims);
bool object_exists(hid_t object_id, const char* name);
void read_dataset(hid_t obj_id, const char* name, hid_t mem_type_id, void* buffer, bool indep);
hid_t open_dataset(hid_t group_id, const char* name);
hid_t open_group(hid_t group_id, const char* name);
}
#endif
