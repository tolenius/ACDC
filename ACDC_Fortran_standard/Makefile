FC = gfortran
FCFLAGS = -O3 -fcheck=bounds -finit-local-zero

run = run_acdc_J_example.f90
get_J = get_acdc_J.f90
driver = driver_acdc_J.f90

eqs = acdc_equations_AN_ions_example.f90
system = acdc_system_AN_ions_example.f90


run: run_acdc_J.o get_acdc_J.o driver.o acdc_equations.o acdc_system.o acdc_simulation_setup.o solution_settings.o vode.o
	$(FC) $(FCFLAGS) $^ -o $@


run_acdc_J.o: $(run) get_acdc_J.o acdc_simulation_setup.o
	$(FC) $(FCFLAGS) -c $< -o $@

get_acdc_J.o: $(get_J) driver.o acdc_system.o acdc_simulation_setup.o
	$(FC) $(FCFLAGS) -c $< -o $@

driver.o: $(driver) acdc_system.o acdc_simulation_setup.o solution_settings.o
	$(FC) $(FCFLAGS) -c $< -o $@

acdc_equations.o: $(eqs) acdc_simulation_setup.o
	$(FC) $(FCFLAGS) -c $< -o $@

acdc_system.o: $(system)
	$(FC) $(FCFLAGS) -c $< -o $@

acdc_simulation_setup.o: acdc_simulation_setup.f90 acdc_system.o
	$(FC) $(FCFLAGS) -c $< -o $@

solution_settings.o: solvers/solution_settings.f90
	$(FC) $(FCFLAGS) -c $< -o $@
    
vode.o: solvers/dvode.f
	$(FC) $(FCFLAGS) -std=legacy -c $< -o $@

.PHONY: clean

clean:
	rm -f *.o *.mod run
