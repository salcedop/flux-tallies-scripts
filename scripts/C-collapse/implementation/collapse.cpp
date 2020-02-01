#include "hdf5_interface.h"
#include <array>
#include <cstring>
#include <sstream>
#include <string>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "string.h"

using namespace std;
using namespace hdf5;

//string filename = "O16.h5";
//extern "C" int strlen(const char*);
extern size_t strlen (const char *__s);
//char* concat(string,string);
int main()
{
hid_t group_id;
hid_t main_group;
hid_t xs;
hsize_t size;
char name[] = "U235";
char grp_name[] = {"U235/(n,2n)/groupr"};
const char* name_ptr = &name[0];
//char* grp_name[] = {"O16/(n,2n)/groupr"};
string FILE_NAME = "/home/salcedop/lib-piecewise-final-practive/MG-lib-save/U235/U235.h5";

//const H5std_string FILE_NAME(path);

hid_t file_id = hdf5::file_open(FILE_NAME,'r');
//hdf5::get_groups(file_id,&grp_name);
group_id = hdf5::open_group(file_id,name_ptr);
hdf5::get_name(group_id,name);
main_group = hdf5::open_group(group_id,"reactions/(n,gamma)");
xs = hdf5::open_dataset(main_group,"groupr");
hdf5::get_shape(xs,&size);
float *xs_data = new float[size];
bool idep = false;
hid_t mem_type = H5T_NATIVE_FLOAT;
hdf5::read_dataset(main_group,"groupr",mem_type,&xs_data[0],idep);

for (hid_t i = 0; i<size; i++){

printf("%f\n",*(xs_data+i));

}

//hid_t group = hdf5::open_group(main_group,grp_name[0]);
return 0;
}
/*
char* concat(string s1, string s2)
{
char* result = new(strlen(s1) + strlen(s2) + 1);
        strcpy(result, s1);
        strcat(result, s2);
        return result;
}
*/
