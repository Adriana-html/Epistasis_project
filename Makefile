#compiladores
CXX = g++
MPICXX = mpic++

# Banderas de compilación
# -O3: Optimización máxima de velocidad
CXXFLAGS = -O3 

#Nombres de los archivos fuente y ejecutables
serial_executable = serial_epistasis
mpi_executable = mpi_epistasis

#ruta de los archivos fuente
serial_src = src/serial_epistasis_v2.cpp
MPI_src = src/mpi_epistasis.cpp

#Reglas de compilación
all: $(serial_executable) $(mpi_executable)

#compilación de la versión serial
$(serial_executable): $(serial_src)
	@echo "Compilando versión serial..."
	$(CXX) $(CXXFLAGS) $(serial_src) -o $(serial_executable)
	@echo "Compilación de la versión serial completada."

#compilación de la versión MPI
$(mpi_executable): $(MPI_src)
	@echo "Compilando versión MPI..."
	$(MPICXX) $(CXXFLAGS) $(MPI_src) -o $(mpi_executable)
	@echo "Compilación de la versión MPI completada."	

#Regla para limpiar los archivos ejecutables
clean:
	@echo "Limpiando archivos ejecutables..."
	rm -f $(serial_executable) $(mpi_executable)
	@echo "Limpieza completada."
	rm -f results/*.txt
	rm -f results/*.csv
	@echo "Limpieza de archivos de resultados completada."

#regla de ayuda para recordar los comandos disponibles
help:
	@echo "Comandos disponibles:"
	@echo "  make        - Compila ambas versiones (serial y MPI)"
	@echo "  make clean  - Elimina los archivos ejecutables y de resultados"
	@echo "  make help   - Muestra esta ayuda"