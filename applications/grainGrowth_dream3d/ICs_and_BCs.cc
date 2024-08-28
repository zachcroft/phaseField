// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::setInitialCondition(const dealii::Point<dim> &p,
                                            const unsigned int        index,
                                            double                   &scalar_IC,
                                            dealii::Vector<double>   &vector_IC,
                                            const std::vector<double> &data)
{
  // ---------------------------------------------------------------------
  // ENTER THE INITIAL CONDITIONS HERE
  // ---------------------------------------------------------------------
  // Enter the function describing conditions for the fields at point "p".
  // Use "if" statements to set the initial condition for each variable
  // according to its variable index

  if (index == 0)
  {
    int Nx = userInputs.domain_size[0];
    int Ny = userInputs.domain_size[1];
    int Nz = userInputs.domain_size[2];

    // Convert cartesion point into flattened array index
    int index = int(p[0] + (Ny+1) * p[1] + (Ny+1)*(Nz+1) * p[2]);

    // Extract the grain ID from the data array
    double data_value = data[index];

    //double data_value = 0.0;
    std::cout << "Data value at index " << index << " is " << data_value << std::endl;
    scalar_IC = data_value;
  }



  // --------------------------------------------------------------------------
}

// ===========================================================================
// FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
// ===========================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::setNonUniformDirichletBCs(const dealii::Point<dim> &p,
                                                  const unsigned int        index,
                                                  const unsigned int        direction,
                                                  const double              time,
                                                  double                   &scalar_BC,
                                                  dealii::Vector<double>   &vector_BC)
{
  // --------------------------------------------------------------------------
  // ENTER THE NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS HERE
  // --------------------------------------------------------------------------
  // Enter the function describing conditions for the fields at point "p".
  // Use "if" statements to set the boundary condition for each variable
  // according to its variable index. This function can be left blank if there
  // are no non-uniform Dirichlet boundary conditions. For BCs that change in
  // time, you can access the current time through the variable "time". The
  // boundary index can be accessed via the variable "direction", which starts
  // at zero and uses the same order as the BC specification in parameters.in
  // (i.e. left = 0, right = 1, bottom = 2, top = 3, front = 4, back = 5).

  // -------------------------------------------------------------------------
}
