cd build
rm -r *
cmake .. -Dlibcloudph++_DIR=/home/pzmij/biblioteki/local_folder/Coal/share/libcloudph++ -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=/home/pzmij/biblioteki/ -DCMAKE_CXX_COMPILER=c++ 
make
