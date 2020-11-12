cd build
rm -r *
##cmake .. -Dlibcloudph++_DIR=/mnt/local/pzmij/UWLCM_singu/singularity/local_folder/multi/share/libcloudph++
echo `pwd`
cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=/mnt/local/pzmij/UWLCM_singu/singularity/local_folder -DCMAKE_CXX_COMPILER=c++ -Dlibcloudph++_DIR=/mnt/local/pzmij/UWLCM_singu/singularity/local_folder/multi/share/libcloudph++
VERBOSE=1 make

