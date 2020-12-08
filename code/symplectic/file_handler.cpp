#include "file_handler.h"

namespace fh{

/**
 * Load data from a file into a c++ standar library vector of type double. Empty values in file are 
 * filled with (quiet) nan.
 * 
 * @param filename: file to be loaded.
 * @param data: vector to be loaded with data from the file.
 * @param delim: file data delimiter. Default: ','.
 * @param  comment: file comment indicator. Default: '#'.
 * 
 * @raises FileNotOpened
 */
void vec_load_file(std::string filename, vector2D(double) &data, char delim, char comment){
    std::ifstream file(filename);
    data.clear();

    if(file.is_open()){
        std::string line;
        while(getline(file, line)){
            std::string temp;
            vector(double) row;
            bool skipped = false;

            for(unsigned int i=0; i<line.size(); i++){
                char current = line[i];
                if (current == comment){
                    skipped = true;
                    break;
                }
                else if(current != delim) temp += line[i];
                else{
                    row.push_back(std::stod(temp));
                    temp = "";
                }
            }

            if(skipped == false){
                row.push_back(std::stod(temp));
                data.push_back(row);
            }
        }
    }
    else throw FileNotOpened();

    file.close();

    __check_vec(data);
}

/**
 * Load data from a file into an array of doubles. Empty values in file are filled with (quiet) nan.
 * Memory management over param 'data' is the user's responsability (However, a function is provided 
 * in this module to aid with this. See fh::clear.)
 * 
 * @param filename: file to be loaded.
 * @param data: pointer to be loaded with data from the file.
 * @param rows: integer where the array's number of rows will be stored.
 * @param columns: integer where the array's number of columns will be stored.
 * @param delim: file data delimiter. Default: ','.
 * @param  comment: file comment indicator. Default: '#'.
 * 
 * @raises FileNotOpened (via vec_load_file)
 */
void load_file(std::string filename, double **&data, int &rows, int &columns, char delim, char comment){
    vector2D(double) data_vec;
    vec_load_file(filename, data_vec, delim, comment);

    int n_rows = data_vec.size(), n_columns = data_vec[0].size();

    data = new double *[n_rows];
    for(int i=0; i<n_rows; i++)
        data[i] = new double[n_columns];

    for(int i=0; i<n_rows; i++)
        for(int j=0; j<n_columns; j++){
            data[i][j] = data_vec[i][j];
        }

    rows = n_rows; columns = n_columns;

}

/**
 * Simplify memory management, freeing memory and reseting pointer to NULL.
 * 
 * @param data: pointer to 2D array to be freed.
 * @param rows: number of rows of the array.
 * @param set_null: define if param 'data' is to be set to null (In case you wish to reuse it).
 */
void clear(double **&data, int rows, bool set_null){
    for(int i=0; i<rows; i++) delete[] data[i];
    delete[] data;

    if(set_null == true) data = NULL;
}

/**
 * Checks that a given vector has every row with the same size. If not, adds append to it.
 */
void __check_vec(vector2D(double) &data){
    int row_len = data[0].size();

    for(unsigned int i=1; i<data.size(); i++){
        int current = data[i].size();

        if(current == row_len) continue;
        else if(current < row_len){
            while(current < row_len){
                data[i].push_back(std::nan("1"));
                current += 1;
            }
            std::cerr << "Missing value in file, NaN added to data. Line " 
            << i+1 << " too short" << std::endl;
        }
        else{
            for(unsigned int j=0; j<i; j++){
                int temp = row_len;
                while(temp < current){
                    data[j].push_back(std::nan("1"));
                    temp += 1;
                }
            }
            row_len = current;
            std::cerr << "Missing values in file, NaN added to data. Line " 
            << i+1 << " too long." << std::endl;
        }
    }
}

/* Get full file as a string */
std::string get_full_file(std::ifstream f){
    std::stringstream buffer;
    buffer << f.rdbuf();
    return buffer.str();
}

} // namespace fh

#undef vector2D
#undef vector
