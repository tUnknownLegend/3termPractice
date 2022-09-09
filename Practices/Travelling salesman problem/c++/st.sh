rm -r ./build/*
cd build
scan-build cmake ..
scan-build make

printf "\n--------------\nSLOW:\n"
./1_rec
printf "\n--------------\nQUICK:\n"
./1_tsp

#valgrind --tool=memcheck --leak-check=yes --error-exitcode=1 ./1_rec
#valgrind --tool=memcheck --leak-check=yes --error-exitcode=1 ./1_tsp

#printf "main:\n--------------\n"
#valgrind --tool=memcheck --leak-check=yes --error-exitcode=1 ./1
#printf "\n--------------\nTests:\n"
#valgrind --tool=memcheck --leak-check=yes --error-exitcode=1 ./Tests