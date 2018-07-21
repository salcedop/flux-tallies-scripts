#include "hdf5_interface.h"

//#include "error.h"
#include <array>
#include <cstring>
#include <sstream>
#include <string>
#include <iostream>
#include "hdf5.h"
#include "hdf5_hl.h"

using namespace std;
namespace hdf5 {

bool
attribute_exists(hid_t obj_id, const char* name)
{
  htri_t out = H5Aexists_by_name(obj_id, ".", name, H5P_DEFAULT);
  return out > 0;
}

bool
object_exists(hid_t object_id, const char* name)
{
  htri_t out = H5LTpath_valid(object_id, name, true);
  if (out < 0) {
    stringstream err_msg;
    cout << "Failed to check if object \"" << name << "\" exists.";
    //error::fatal_error(err_msg);
    //cout << err_msg;
  }
  return (out > 0);
}

hid_t
open_dataset(hid_t group_id, const char* name)
{
  if (object_exists(group_id, name)) {
    return H5Dopen(group_id, name, H5P_DEFAULT);
  } else {
    stringstream err_msg;
    cout << "Group \"" << name << "\" does not exist";
    //error::fatal_error(err_msg);
    //cout << err_msg;
  }
}

hid_t
file_open(const char* filename, char mode)
{
  bool create;
  unsigned int flags;
  switch (mode) {
    case 'r':
    case 'a':
      create = false;
      flags = (mode == 'r' ? H5F_ACC_RDONLY : H5F_ACC_RDWR);
      break;
    case 'w':
    case 'x':
      create = true;
      flags = (mode == 'x' ? H5F_ACC_EXCL : H5F_ACC_TRUNC);
      break;
    default:
      stringstream err_msg;
      cout <<  "Invalid file mode: " << mode;
      //error::fatal_error(err_msg);
      //cout << err_msg;
}
  hid_t plist = H5P_DEFAULT;
  hid_t file_id;
  if (create) {
    file_id = H5Fcreate(filename, flags, H5P_DEFAULT,plist);
  } else {
    file_id = H5Fopen(filename, flags,plist);
  }
  if (file_id < 0) {
    stringstream msg;
    msg << "Failed to open HDF5 file with mode '" << mode << "': " << filename;
    //error::fatal_error(msg);
    //cout << msg;
  }

  return file_id;
}

hid_t
file_open(const string& filename, char mode)
{
  file_open(filename.c_str(), mode);
}

void file_close(hid_t file_id)
{
  H5Fclose(file_id);
}


void
get_groups(hid_t group_id, char* name[])
{
  // Determine number of links in the group
  H5G_info_t info;
  H5Gget_info(group_id, &info);

  // Iterate over links to get names
  H5O_info_t oinfo;
  hsize_t count = 0;
  size_t size;
  for (hsize_t i = 0; i < info.nlinks; ++i) {
    // Determine type of object (and skip non-group)
    H5Oget_info_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC, i, &oinfo,
                       H5P_DEFAULT);
    if (oinfo.type != H5O_TYPE_GROUP) continue;

    // Get size of name
    size = 1 + H5Lget_name_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC,
                                  i, nullptr, 0, H5P_DEFAULT);

    // Read name
    H5Lget_name_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_INC, i,
                       name[count], size, H5P_DEFAULT);
    count += 1;
  }
}

hid_t
open_group(hid_t group_id, const char* name)
{
  if (object_exists(group_id, name)) {
    return H5Gopen(group_id, name, H5P_DEFAULT);
  } else {
    stringstream err_msg;
    cout << "Group \"" << name << "\" does not exist";
    //error::fatal_error(err_msg);
    //cout << err_msg;
  }
}

void
get_name(hid_t obj_id, char* name)
{
  size_t size = 1 + H5Iget_name(obj_id, nullptr, 0);
  H5Iget_name(obj_id, name, size);
}

void
get_shape(hid_t obj_id, hsize_t* dims)
{
  auto type = H5Iget_type(obj_id);
  hid_t dspace;
  if (type == H5I_DATASET) {
    dspace = H5Dget_space(obj_id);
  } else if (type == H5I_ATTR) {
    dspace = H5Aget_space(obj_id);
  }
  H5Sget_simple_extent_dims(dspace, dims, nullptr);
  H5Sclose(dspace);
}


void
get_shape_attr(hid_t obj_id, const char* name, hsize_t* dims)
{
  hid_t attr = H5Aopen(obj_id, name, H5P_DEFAULT);
  hid_t dspace = H5Aget_space(attr);
  H5Sget_simple_extent_dims(dspace, dims, nullptr);
  H5Sclose(dspace);
  H5Aclose(attr);
}

void
read_dataset(hid_t obj_id, const char* name, hid_t mem_type_id,
             void* buffer, bool indep)
{
  hid_t dset = obj_id;
  if (name) dset = open_dataset(obj_id, name);

    H5Dread(dset, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
  
  if (name) H5Dclose(dset);
}
}
