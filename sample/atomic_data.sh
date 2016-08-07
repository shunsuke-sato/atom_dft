# Hydrogen
cd H
../../atom_dft < input_H > log_H.log
gnuplot plot.plt
cd ..

# Helium
cd He
../../atom_dft < input_He > log_He.log
gnuplot plot.plt
cd ..

# Lithium
cd Li
../../atom_dft < input_Li > log_Li.log
gnuplot plot.plt
cd ..

# Berylium
cd Be
../../atom_dft < input_Be > log_Be.log
gnuplot plot.plt
cd ..

# Carbon
cd C
../../atom_dft < input_C > log_C.log
gnuplot plot.plt
cd ..

# Nitrogen
cd N
../../atom_dft < input_N > log_N.log
gnuplot plot.plt
cd ..

# Oxygen
cd O
../../atom_dft < input_O > log_O.log
gnuplot plot.plt
cd ..

# Fluorine
cd F
../../atom_dft < input_F > log_F.log
gnuplot plot.plt
cd ..

# Neon
cd Ne
../../atom_dft < input_Ne > log_Ne.log
gnuplot plot.plt
cd ..
