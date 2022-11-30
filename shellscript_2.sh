
ifort -o ami.x LangDyn_2D_Inter_DD_1k.f90
./ami.x
ifort -O3 -o ami.x LangDyn_2D_Inter_DD_1k.f90
./ami.x
ifort -o ami.x LangDyn_2D_Inter_DD_omp_1k.f90 -fopenmp
./ami.x
ifort -O3 -o ami.x LangDyn_2D_Inter_DD_omp_1k.f90 -fopenmp
./ami.x
export OMP_NUM_THREADS=2
./ami.x
export OMP_NUM_THREADS=3
./ami.x
export OMP_NUM_THREADS=4
./ami.x

ifort -o ami.x LangDyn_2D_Inter_DD_2k.f90
./ami.x
ifort -O3 -o ami.x LangDyn_2D_Inter_DD_2k.f90
./ami.x
ifort -o ami.x LangDyn_2D_Inter_DD_omp_2k.f90 -fopenmp
./ami.x
ifort -O3 -o ami.x LangDyn_2D_Inter_DD_omp_2k.f90 -fopenmp
./ami.x
export OMP_NUM_THREADS=2
./ami.x
export OMP_NUM_THREADS=3
./ami.x
export OMP_NUM_THREADS=4
./ami.x

ifort -o ami.x LangDyn_2D_Inter_DD_3k.f90
./ami.x
ifort -O3 -o ami.x LangDyn_2D_Inter_DD_3k.f90
./ami.x
ifort -o ami.x LangDyn_2D_Inter_DD_omp_3k.f90 -fopenmp
./ami.x
ifort -O3 -o ami.x LangDyn_2D_Inter_DD_omp_3k.f90 -fopenmp
./ami.x
export OMP_NUM_THREADS=2
./ami.x
export OMP_NUM_THREADS=3
./ami.x
export OMP_NUM_THREADS=4
./ami.x

ifort -o ami.x LangDyn_2D_Inter_DD_4k.f90
./ami.x
ifort -O3 -o ami.x LangDyn_2D_Inter_DD_4k.f90
./ami.x
ifort -o ami.x LangDyn_2D_Inter_DD_omp_4k.f90 -fopenmp
./ami.x
ifort -O3 -o ami.x LangDyn_2D_Inter_DD_omp_4k.f90 -fopenmp
./ami.x
export OMP_NUM_THREADS=2
./ami.x
export OMP_NUM_THREADS=3
./ami.x
export OMP_NUM_THREADS=4
./ami.x
