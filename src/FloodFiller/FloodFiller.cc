#include "../../include/FloodFiller.h"

#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <vector>

// Function to calculate sets of grains described by
// a single order parameter
template <int dim, int degree>
void
FloodFiller<dim, degree>::calcGrainSets_optimal(dealii::FESystem<dim>   &fe,
                                                dealii::DoFHandler<dim> &dof_handler,
                                                vectorType              *solution_field,
                                                double                   threshold_lower,
                                                double                   threshold_upper,
                                                unsigned int order_parameter_index,
                                                std::vector<GrainSet<dim>> &grain_sets)
{
  // Print statement for debugging
  // pcout << "Entering calcGrainSets function call..." << std::endl;

  // Create variable di, which is the current cell
  typename dealii::DoFHandler<dim>::cell_iterator di = dof_handler.begin(level);

  // pcout << "Creating grain_set and pushing back..." << std::endl;
  // GrainSet<dim> grain_set;
  // grain_sets.push_back(grain_set);
  // grain_sets.back().setOrderParameterIndex(order_parameter_index);
  // grain_sets.back().setGrainIndex(22);

  unsigned int current_grain_ID;
  unsigned int num_grains_found = 0;

  std::vector<unsigned int> grain_IDs_found;

  // Loop over all cells, locating grains
  unsigned int numberOfCellsIterated = 0;
  while (di != dof_handler.end(level))
    {
      // Print the cell iteration
      // pcout << "di = " << di << "\n";

      if (di->has_children())
        {
          pcout << "current cell has children\n";
        }
      else
        {
          // Get the grain ID from the current cell if it hasn't been checked already
          current_grain_ID = checkCell_optimal(di,
                                               dof_handler.end(level),
                                               solution_field,
                                               threshold_lower,
                                               threshold_upper);

          // Only process the data if the cell has an ID greater than the minimum allowed
          // (e.g., zero)
          // if (current_grain_ID > 0)
          // if (current_grain_ID == 30)
          if (current_grain_ID > 0)
            {
              std::cout << "Current GID: " << current_grain_ID << "\n";

              // Set a threshold range for flood filling the current grain
              threshold_lower = current_grain_ID - 0.01;
              threshold_upper = current_grain_ID + 0.01;

              // Create a new grain set object to fill
              GrainSet<dim> new_grain_set;
              new_grain_set.setGrainIndex(current_grain_ID);
              new_grain_set.setOrderParameterIndex(order_parameter_index);
              grain_sets.push_back(new_grain_set);

              // Fill in the grain
              bool grain_assigned = false;
              queueFloodFill<typename dealii::DoFHandler<dim>::cell_iterator>(
                di,
                dof_handler.end(level),
                solution_field,
                threshold_lower,
                threshold_upper,
                grain_sets,
                grain_assigned);

              // pcout << "grain_assigned = " << grain_assigned << "\n";

              num_grains_found++;

              /*
              if (grain_assigned == true && num_grains_found < 1)
                {
                  grain_IDs_found.push_back(current_grain_ID);
                  num_grains_found++;

                  // GrainSet<dim> new_grain_set;
                  // new_grain_set.setOrderParameterIndex(order_parameter_index);
                  //  Call the GrainSet method to label the grain ID
                  //  new_grain_set.setGrainIndex(current_grain_ID);

                  // Push back this grain
                  // grain_sets.push_back(new_grain_set);
                  // break; // Zach added to see which grain is first
                }
                */
            }
        }

      ++di;
      ++numberOfCellsIterated;
    }

  // Generate global list of the grains, merging grains split between multiple
  // processors

  if (dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) > 1)
    {
      // Send the grain set info to all processors so everyone has the full list
      // createGlobalGrainSetList(grain_sets);

      // Merge grains that are split across processors
      // mergeSplitGrains(grain_sets);
    }

  // Check the grain IDs that are stored
  pcout << "Check grain IDs that are stored in grain_sets"
        << "\n";
  for (unsigned int g = 0; g < grain_sets.size() + 1; g++)
    {
      std::cout << "g = " << g << "\n";
      std::cout << "grain_sets[g].getGrainIndex(): " << grain_sets[g].getGrainIndex()
                << "\n";
    }

  // Print out the list of grains found & filled
  std::cout << "Number of grains found: " << num_grains_found << "\n";
  std::cout << "List of grains found:"
            << "\n";
  for (unsigned int GID : grain_IDs_found)
    {
      std::cout << GID << " ";
    }
  std::cout << "\n";

  std::cout << "numberOfCellsIterated = " << numberOfCellsIterated << "\n";
}

//
//
//
// Function to calculate sets of grains described by
// a single order parameter
template <int dim, int degree>
void
FloodFiller<dim, degree>::calcGrainSets(dealii::FESystem<dim>      &fe,
                                        dealii::DoFHandler<dim>    &dof_handler,
                                        vectorType                 *solution_field,
                                        double                      threshold_lower,
                                        double                      threshold_upper,
                                        unsigned int                order_parameter_index,
                                        std::vector<GrainSet<dim>> &grain_sets)
{
  // Print statement for debugging
  // pcout << "Entering calcGrainSets function call..." << std::endl;

  // Create variable di, which is the current cell
  typename dealii::DoFHandler<dim>::cell_iterator di = dof_handler.begin(level);

  pcout << "dof_handler.end(level) = " << dof_handler.end(level) << "\n";

  // Loop through the whole mesh and set the user flags to false
  // (so everything is considered unmarked)
  // pcout << "Looping through the mesh to clear flags..." << std::endl;
  /*
  while (di != dof_handler.end(level))
    {
      if (!di->has_children())
        {
          // Clear each cell of the 'marked' marking
          di->clear_user_flag();
        }
      ++di;
    }
  */
  // Zach: commented out the clearing of flags for visited cells (seemed to have no
  // effect)

  // Print statement for debugging
  // std::cout << "Finished clearing cell markers..." << std::endl;

  // pcout << "Creating grain_set and pushing back..." << std::endl;
  GrainSet<dim> grain_set;
  grain_sets.push_back(grain_set);
  grain_sets.back().setOrderParameterIndex(order_parameter_index);

  // The flood fill loop
  di                                 = dof_handler.begin(level);
  unsigned int numberOfCellsIterated = 0;

  // bool grain_assigned = false; // Zach added this line for testing (might break the
  // code)

  // Note: 'level' refers to refine factor
  while (di != dof_handler.end(level)) // && !grain_assigned)
    {
      // Print statement for debugging
      // pcout << "Top of while loop iteration over cells...\n";
      // if (di->has_children())
      //  {
      //    pcout << "current cell has children\n";
      //  }

      if (!di->has_children())
        {
          // pcout << "Inside if(!di->has_children) section...\n";
          // pcout << "numberOfCellsIterated = " << numberOfCellsIterated << "\n";

          bool grain_assigned = false;

          // Fill in the current grain
          queueFloodFill<typename dealii::DoFHandler<dim>::cell_iterator>(di,
                                                                          dof_handler.end(
                                                                            level),
                                                                          solution_field,
                                                                          threshold_lower,
                                                                          threshold_upper,
                                                                          grain_sets,
                                                                          grain_assigned);

          // Check if the grain has been assigned to a grain set.
          // If it hasn't, add it to a new grain set.
          if (grain_assigned)
            {
              GrainSet<dim> new_grain_set;
              new_grain_set.setOrderParameterIndex(order_parameter_index);
              grain_sets.push_back(new_grain_set);
            }
        }

      ++di;
      ++numberOfCellsIterated;
    }

  // Generate global list of the grains, merging grains split between multiple
  // processors
  if (dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) > 1)
    {
      // Send the grain set info to all processors so everyone has the full list
      createGlobalGrainSetList(grain_sets);

      // Merge grains that are split across processors
      mergeSplitGrains(grain_sets);
    }
}

// Function for performing the flood fill operation
template <int dim, int degree>
template <typename T>
void
FloodFiller<dim, degree>::queueFloodFill(T                           di,
                                         T                           di_end,
                                         vectorType                 *solution_field,
                                         double                      threshold_lower,
                                         double                      threshold_upper,
                                         std::vector<GrainSet<dim>> &grain_sets,
                                         bool                       &grain_assigned)
{
  pcout << "Entering queueFloodFill()..." << std::endl;
  pcout << "threshold_lower = " << threshold_lower << "\n";
  pcout << "threshold_upper = " << threshold_upper << "\n";

  // Check if checkCell is false (i.e., if the cell has already been visited or the value
  // isn't the grain ID is not the one we're currently trying to Locate then skip this
  // iteration of queueFloodFill)
  if (!checkCell(di, di_end, solution_field, threshold_lower, threshold_upper))
    {
      // pcout << "!checkCell is true...\n";
      return;
    }

  // If the code makes it to this point, then checkCell has returned true and we found a
  // starting point for current the grain we're trying to locate. First, we will mark this
  // cell as visited by setting a flag.
  di->set_user_flag();

  // Print the starting location of the grain for flood fill
  pcout << "Grain of interest found at: " << di->vertex(0) << "\n";

  // Make a queue for the flood fill
  std::queue<T> floodQueue;
  floodQueue.emplace(di);

  // Add the vertices of the current cell to a list for this grain since we know it exists
  // at the current cell
  std::vector<dealii::Point<dim>> vertex_list;
  for (unsigned int v = 0; v < dealii::Utilities::fixed_power<dim>(2.0); v++)
    {
      vertex_list.push_back(di->vertex(v));
    }
  grain_sets.back().addVertexList(vertex_list);
  vertex_list.clear();

  // Counting variable for # of cells visited for the current filling operation
  unsigned int numberOfCellsInFill = 1;

  // Queue Loop for flood fill:
  while (floodQueue.size() > 0)
    {
      // De-queue the current cell
      T currentCell = floodQueue.front();
      floodQueue.pop();

      // Loop over the current cell's neighbors since they are possible candidates for the
      // current grain
      std::vector<T> possible_neighbors;
      for (unsigned int n = 0; n < 2 * dim; n++)
        {
          // Check if the neighbor cell contains the desired value (i.e., it is connected
          // to the grain we care about). If it is, add it to the possible_neighbors
          // vector
          if (checkCell(currentCell->neighbor(n),
                        di_end,
                        solution_field,
                        threshold_lower,
                        threshold_upper))
            {
              possible_neighbors.push_back(currentCell->neighbor(n));
            }
        }

      // Loop over the neighbor list that contained the desired grain ID
      for (unsigned int i = 0; i < possible_neighbors.size(); i++)
        {
          // Flag these neighbor cells as visited and add their vertices to the index list
          // for the grain set
          possible_neighbors[i]->set_user_flag();
          std::vector<dealii::Point<dim>> vertex_list;
          for (unsigned int v = 0; v < dealii::Utilities::fixed_power<dim>(2.0); v++)
            {
              vertex_list.push_back(possible_neighbors[i]->vertex(v));
            }
          grain_sets.back().addVertexList(vertex_list);

          // Enqueue (i.e., add the neighbor cells to the floodQueue list so that their
          // neighbors can be checked recursively)
          floodQueue.emplace(possible_neighbors[i]);

          ++numberOfCellsInFill;
        }
    }

  pcout << "Grain filled: found " << numberOfCellsInFill << " cells in grain.\n";
  // Now that the flood fill has been completed for the current grain, flag this grain as
  // being assigned
  grain_assigned = true;
}

template <int dim, int degree>
template <typename T>
bool
FloodFiller<dim, degree>::checkCell(T           di,
                                    T           di_end,
                                    vectorType *solution_field,
                                    double      threshold_lower,
                                    double      threshold_upper)
{
  // Print statement for debugging
  // std::cout << "Entering checkCell function call..." << std::endl;

  if (di != di_end && !di->user_flag_set() && !di->has_children() &&
      di->is_locally_owned())
    {
      // Print statement for debugging
      // std::cout << "Inside checkCell if statement..." << std::endl;

      std::vector<double> var_values(num_quad_points);
      fe_values.reinit(di);
      fe_values.get_function_values(*solution_field, var_values);

      /*
      // Older version for average value
      double ele_val = 0.0;
      for (unsigned int q_point=0; q_point < num_quad_points; ++q_point){
        for (unsigned int i=0; i < dofs_per_cell; ++i){
          ele_val += fe_values.shape_value(i, q_point) * var_values[q_point] *
      quadrature.weight(q_point);
        }
      }
      std::cout << "ele_val: " << ele_val << std::endl;
      if (threshold_lower < ele_val && ele_val < threshold_upper){
        return true;
      }
      */

      // Newer version for most frequent value
      double                ele_val;
      std::map<double, int> quadratureValues;
      int                   maxNumberSeen         = 0;
      double                mostCommonQPointValue = 0;

      // Loop over quadrature points
      for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
        {
          // Print statement for debugging
          // std::cout << "Top of checkCell for loop iteration..." << std::endl;

          if (var_values[q_point] > 0)
            {
              ++quadratureValues[var_values[q_point]];
            }
          if (quadratureValues[var_values[q_point]] > maxNumberSeen)
            {
              maxNumberSeen         = quadratureValues[var_values[q_point]];
              mostCommonQPointValue = var_values[q_point];
            }
        }

      ele_val = mostCommonQPointValue;
      // pcout << "ele_val: " << ele_val << std::endl;

      // Zach modified to have equals sign
      if (threshold_lower <= ele_val && ele_val <= threshold_upper)
        {
          // pcout << "checkCell is true (!checkCell is false)...\n";
          // pcout << "threshold lower = " << threshold_lower << "\n";
          // pcout << "threshold upper = " << threshold_upper << "\n";
          return true;
        }
    }

  // Print statement for debugging
  // pcout << "checkCell is false (!checkCell is true)...\n";

  return false;
}

template <int dim, int degree>
template <typename T>
double
FloodFiller<dim, degree>::checkCell_optimal(T           di,
                                            T           di_end,
                                            vectorType *solution_field,
                                            double      threshold_lower,
                                            double      threshold_upper)
{
  // Print statement for debugging
  // std::cout << "Entering checkCell function call..." << std::endl;

  if (di != di_end && !di->user_flag_set() && !di->has_children() &&
      di->is_locally_owned())
    {
      // Print statement for debugging
      // std::cout << "Inside checkCell if statement..." << std::endl;

      std::vector<double> var_values(num_quad_points);
      fe_values.reinit(di);
      fe_values.get_function_values(*solution_field, var_values);

      // Newer version for most frequent value
      double                ele_val;
      std::map<double, int> quadratureValues;
      int                   maxNumberSeen         = 0;
      double                mostCommonQPointValue = 0;

      // Loop over quadrature points
      for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
        {
          // Print statement for debugging
          // std::cout << "Top of checkCell for loop iteration..." << std::endl;

          if (var_values[q_point] > 0)
            {
              ++quadratureValues[var_values[q_point]];
            }
          if (quadratureValues[var_values[q_point]] > maxNumberSeen)
            {
              maxNumberSeen         = quadratureValues[var_values[q_point]];
              mostCommonQPointValue = var_values[q_point];
            }
        }

      ele_val = mostCommonQPointValue;
      // pcout << "ele_val: " << ele_val << std::endl;

      // Return the grain ID
      return ele_val;

      // if (threshold_lower < ele_val && ele_val < threshold_upper)
      // {
      //   return true;
      // }
    }

  return false;
}

//
//
//
//
//
//
//
//
//
// Function for performing the flood fill operation.
template <int dim, int degree>
template <typename T>
void
FloodFiller<dim, degree>::recursiveFloodFill(T                           di,
                                             T                           di_end,
                                             vectorType                 *solution_field,
                                             double                      threshold_lower,
                                             double                      threshold_upper,
                                             unsigned int               &grain_index,
                                             std::vector<GrainSet<dim>> &grain_sets,
                                             bool                       &grain_assigned)
{
  if (di != di_end)
    {
      // Check if the cell has been marked yet
      bool cellMarked = di->user_flag_set();

      if (!cellMarked)
        {
          if (di->has_children())
            {
              // Call recursiveFloodFill on the element's children
              for (unsigned int n = 0; n < di->n_children(); n++)
                {
                  recursiveFloodFill<T>(di->child(n),
                                        di_end,
                                        solution_field,
                                        threshold_lower,
                                        threshold_upper,
                                        grain_index,
                                        grain_sets,
                                        grain_assigned);
                }
            }
          else
            {
              if (di->is_locally_owned())
                {
                  di->set_user_flag();

                  dealii::FEValues<dim> fe_values(*fe, quadrature, dealii::update_values);
                  std::vector<double>   var_values(num_quad_points);
                  std::vector<dealii::Point<dim>> q_point_list(num_quad_points);

                  // Get the most common value for the element
                  fe_values.reinit(di);
                  fe_values.get_function_values(*solution_field, var_values);

                  double                ele_val;
                  std::map<double, int> quadratureValues;
                  int                   maxNumberSeen         = 0;
                  double                mostCommonQPointValue = -1;
                  for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
                    {
                      // Add the number of times that var_values[q_point] has been seen
                      ++quadratureValues[var_values[q_point]];
                      if (quadratureValues[var_values[q_point]] > maxNumberSeen)
                        {
                          maxNumberSeen         = quadratureValues[var_values[q_point]];
                          mostCommonQPointValue = var_values[q_point];
                        }
                    }
                  ele_val = mostCommonQPointValue;

                  if (ele_val > threshold_lower && ele_val < threshold_upper)
                    {
                      grain_assigned = true;

                      std::vector<dealii::Point<dim>> vertex_list;
                      for (unsigned int v = 0;
                           v < dealii::Utilities::fixed_power<dim>(2.0);
                           v++)
                        {
                          vertex_list.push_back(di->vertex(v));
                        }
                      grain_sets.back().addVertexList(vertex_list);

                      // Call recursiveFloodFill on the element's neighbors
                      for (unsigned int n = 0; n < 2 * dim; n++)
                        {
                          recursiveFloodFill<T>(di->neighbor(n),
                                                di_end,
                                                solution_field,
                                                threshold_lower,
                                                threshold_upper,
                                                grain_index,
                                                grain_sets,
                                                grain_assigned);
                        }
                    }
                }
            }
        }
    }
}

// =================================================================================
// All-to-all communication of the grain sets
// =================================================================================
template <int dim, int degree>
void
FloodFiller<dim, degree>::createGlobalGrainSetList(
  std::vector<GrainSet<dim>> &grain_sets) const
{
  int numProcs = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  // int thisProc=dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  unsigned int num_grains_local = grain_sets.size();

  // Convert the grain_set object into a group of vectors
  std::vector<unsigned int> order_parameters;
  std::vector<unsigned int> num_elements;
  std::vector<double>       vertices;

  for (unsigned int g = 0; g < grain_sets.size(); g++)
    {
      order_parameters.push_back(grain_sets.at(g).getOrderParameterIndex());

      std::vector<std::vector<dealii::Point<dim>>> vertex_list =
        grain_sets[g].getVertexList();
      num_elements.push_back(vertex_list.size());

      for (unsigned int c = 0; c < num_elements[g]; c++)
        {
          for (unsigned int v = 0; v < dealii::Utilities::fixed_power<dim>(2.0); v++)
            {
              for (unsigned int d = 0; d < dim; d++)
                {
                  vertices.push_back(vertex_list[c][v][d]);
                }
            }
        }
    }

  unsigned int num_vertices = 0;
  for (unsigned int g = 0; g < grain_sets.size(); g++)
    {
      num_vertices += num_elements[g] * dealii::Utilities::fixed_power<dim>(2) * dim;
    }

  // Communicate how many grains each core has
  std::vector<int> num_grains_per_core(numProcs, 0);

  MPI_Allgather(&num_grains_local,
                1,
                MPI_INT,
                &num_grains_per_core[0],
                1,
                MPI_INT,
                MPI_COMM_WORLD);

  int num_grains_global =
    std::accumulate(num_grains_per_core.begin(), num_grains_per_core.end(), 0);

  // Communicate the order_parameters
  std::vector<int> offset(numProcs, 0);
  for (int n = 1; n < numProcs; n++)
    {
      offset[n] = offset[n - 1] + num_grains_per_core[n - 1];
    }

  std::vector<unsigned int> order_parameters_global(num_grains_global, 0);

  MPI_Allgatherv(&order_parameters[0],
                 num_grains_local,
                 MPI_UNSIGNED,
                 &order_parameters_global[0],
                 &num_grains_per_core[0],
                 &offset[0],
                 MPI_UNSIGNED,
                 MPI_COMM_WORLD);

  // Communicate the number of elements
  std::vector<unsigned int> num_elements_global(num_grains_global, 0);

  MPI_Allgatherv(&num_elements[0],
                 num_grains_local,
                 MPI_UNSIGNED,
                 &num_elements_global[0],
                 &num_grains_per_core[0],
                 &offset[0],
                 MPI_UNSIGNED,
                 MPI_COMM_WORLD);

  // Communicate the vertices
  unsigned int total_elements =
    std::accumulate(num_elements_global.begin(), num_elements_global.end(), 0);
  int num_vertices_global = (unsigned int) total_elements *
                            dealii::Utilities::fixed_power<dim>(2) * (unsigned int) dim;
  std::vector<double> vertices_global(num_vertices_global, 0);

  std::vector<int> num_vertices_per_core;

  unsigned int g = 0;
  for (int i = 0; i < numProcs; i++)
    {
      int num_vert_single_core = 0;
      for (int j = 0; j < num_grains_per_core.at(i); j++)
        {
          num_vert_single_core += num_elements_global.at(g) *
                                  dealii::Utilities::fixed_power<dim>(2) *
                                  (unsigned int) dim;
          g++;
        }
      num_vertices_per_core.push_back(num_vert_single_core);
    }

  offset.at(0) = 0;
  for (int n = 1; n < numProcs; n++)
    {
      offset[n] = offset[n - 1] + num_vertices_per_core[n - 1];
    }

  MPI_Allgatherv(&vertices[0],
                 num_vertices,
                 MPI_DOUBLE,
                 &vertices_global[0],
                 &num_vertices_per_core[0],
                 &offset[0],
                 MPI_DOUBLE,
                 MPI_COMM_WORLD);

  // Put the GrainSet objects back together
  grain_sets.clear();

  for (int g = 0; g < num_grains_global; g++)
    {
      GrainSet<dim> new_grain_set;
      for (unsigned int c = 0; c < num_elements_global.at(g); c++)
        {
          std::vector<dealii::Point<dim>> verts;
          for (unsigned int v = 0; v < dealii::Utilities::fixed_power<dim>(2.0); v++)
            {
              double coords[dim];
              for (unsigned int d = 0; d < dim; d++)
                {
                  coords[d] = vertices_global.front();
                  vertices_global.erase(vertices_global.begin());
                }
              dealii::Tensor<1, dim> tensor_coords(coords);
              dealii::Point<dim>     vert(tensor_coords);
              verts.push_back(vert);
            }
          new_grain_set.addVertexList(verts);
        }
      new_grain_set.setOrderParameterIndex(order_parameters_global.at(g));
      grain_sets.push_back(new_grain_set);
    }
}

// =================================================================================
// Check to see if any grains on different processors share vertices
// =================================================================================

template <int dim, int degree>
void
FloodFiller<dim, degree>::mergeSplitGrains(std::vector<GrainSet<dim>> &grain_sets) const
{
  // Loop though each vertex in the base grain "g"
  for (unsigned int g = 0; g < grain_sets.size(); g++)
    {
      std::vector<std::vector<dealii::Point<dim>>> vertex_list =
        grain_sets[g].getVertexList();

      // Now cycle through the other grains to find overlapping elements
      for (unsigned int g_other = g + 1; g_other < grain_sets.size(); g_other++)
        {
          bool matching_vert = false;

          std::vector<std::vector<dealii::Point<dim>>> vertex_list_other =
            grain_sets[g_other].getVertexList();

          for (unsigned int c = 0; c < vertex_list.size(); c++)
            {
              for (unsigned int v = 0; v < dealii::Utilities::fixed_power<dim>(2.0); v++)
                {
                  for (unsigned int c_other = 0; c_other < vertex_list_other.size();
                       c_other++)
                    {
                      for (unsigned int v_other = 0;
                           v_other < dealii::Utilities::fixed_power<dim>(2.0);
                           v_other++)
                        {
                          // Check if the vertices match
                          if (vertex_list[c][v] == vertex_list_other[c_other][v_other])
                            {
                              matching_vert = true;
                              break;
                            }
                          if (matching_vert)
                            {
                              break;
                            }
                        }
                      if (matching_vert)
                        {
                          break;
                        }
                    }
                  if (matching_vert)
                    {
                      break;
                    }
                }
              if (matching_vert)
                {
                  break;
                }
            }

          if (matching_vert)
            {
              for (unsigned int c_base = 0; c_base < vertex_list.size(); c_base++)
                {
                  grain_sets[g_other].addVertexList(vertex_list.at(c_base));
                }
              grain_sets.erase(grain_sets.begin() + g);
              g--;
              break;
            }
        }
    }
}

// Template instantiations
template class FloodFiller<2, 1>;
template class FloodFiller<2, 2>;
template class FloodFiller<2, 3>;
template class FloodFiller<2, 4>;
template class FloodFiller<3, 1>;
template class FloodFiller<3, 2>;
template class FloodFiller<3, 3>;
template class FloodFiller<3, 4>;
