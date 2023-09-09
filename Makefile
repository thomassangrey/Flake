
SRC = ./SRC
OBJ = ./OBJ
BIN = ./BIN

FIpy: flake.o flake_albedo_ref.o flake_configure.o \
src_flake_interface_1D.o flake_paramoptic_ref.o flake_parameters.o \
SfcFlx.o flake_derivedtypes.o data_parameters.o
	gfortran -I$(OBJ) -shared $(SRC)/flake_interface_wrapper.f90 \
	$(SRC)/f2py3_SfcFlx_lwradatm.f90 $(OBJ)/flake.o $(OBJ)/flake_albedo_ref.o \
	$(OBJ)/flake_configure.o $(OBJ)/src_flake_interface_1D.o \
	$(OBJ)/flake_paramoptic_ref.o $(OBJ)/flake_parameters.o $(OBJ)/SfcFlx.o \
	$(OBJ)/flake_derivedtypes.o $(OBJ)/data_parameters.o -o $(BIN)/FI_shared.so

flake_parameters.o: $(SRC)/flake_parameters.f90 data_parameters.o flake_derivedtypes.o
	gfortran -J$(OBJ) -I$(OBJ) -c $(SRC)/flake_parameters.f90 -o $(OBJ)/flake_parameters.o

data_parameters.o: $(SRC)/data_parameters.f90
	gfortran -J$(OBJ) -c $(SRC)/data_parameters.f90 -o $(OBJ)/data_parameters.o

flake_derivedtypes.o: $(SRC)/flake_derivedtypes.f90 data_parameters.o
	gfortran -J$(OBJ) -I$(OBJ) -c $(SRC)/flake_derivedtypes.f90 -o $(OBJ)/flake_derivedtypes.o

flake.o: $(SRC)/flake.f90 data_parameters.o flake_derivedtypes.o \
flake_parameters.o flake_configure.o $(SRC)/flake_driver.incf \
$(SRC)/flake_buoypar.incf $(SRC)/flake_buoypar.incf \
$(SRC)/flake_snowdensity.incf $(SRC)/flake_snowheatconduct.incf
	gfortran -J$(OBJ) -I$(OBJ) -c $(SRC)/flake.f90 -o $(OBJ)/flake.o

flake_albedo_ref.o: $(SRC)/flake_albedo_ref.f90 data_parameters.o
	gfortran -J$(OBJ) -I$(OBJ)  -c $(SRC)/flake_albedo_ref.f90 -o $(OBJ)/flake_albedo_ref.o

flake_configure.o: $(SRC)/flake_configure.f90 data_parameters.o
	gfortran -J$(OBJ) -I$(OBJ) -c $(SRC)/flake_configure.f90 -o $(OBJ)/flake_configure.o

src_flake_interface_1D.o: $(SRC)/src_flake_interface_1D.f90 data_parameters.o \
flake.o flake_parameters.o flake_paramoptic_ref.o SfcFlx.o flake_albedo_ref.o \
flake_derivedtypes.o
	gfortran -I$(OBJ) -c $(SRC)/src_flake_interface_1D.f90 \
	-o $(OBJ)/src_flake_interface_1D.o

flake_paramoptic_ref.o: $(SRC)/flake_paramoptic_ref.f90 data_parameters.o 
	gfortran -J$(OBJ) -I$(OBJ) -c $(SRC)/flake_paramoptic_ref.f90 \
	-o $(OBJ)/flake_paramoptic_ref.o

SfcFlx.o: $(SRC)/SfcFlx.f90 data_parameters.o flake_parameters.o \
$(SRC)/SfcFlx_lwradwsfc.incf $(SRC)/SfcFlx_momsenlat.incf $(SRC)/SfcFlx_lwradatm.incf \
$(SRC)/SfcFlx_rhoair.incf $(SRC)/SfcFlx_momsenlat.incf
	gfortran -J$(OBJ) -I$(OBJ) -c $(SRC)/SfcFlx.f90 -o $(OBJ)/SfcFlx.o