#ifndef Body_HH
#define Body_HH

#include "./Mesh.hh"
#include "./PField.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

namespace PRISMS
{
    /// A class for a Body: a combination of Mesh and Field(s))
    ///
    template< class Coordinate, int DIM>
    class Body
    {
    public:
        Mesh<Coordinate, DIM> mesh;
        std::vector< PField<Coordinate, double, DIM> > scalar_field;
        //std::vector< PField<Coordinate, std::vector<double>, DIM > > vector_field;
        //std::vector< PField<Coordinate, Tensor<double>, DIM > > tensor_field;
        // ----------------------------------------------------------
        // Constructors
        Body(){};

        /// Function for reading from a VTK file
        ///   For now:
        ///      only ASCII files
        ///      only rectilinear grids (though output as UNSTRUCTURED_GRID)
        ///      only (2d) Quad elements
        ///
        void read_vtk( const std::string &vtkfile)
        {
            std::cout << "Begin reading vtk file" << std::endl;
            std::ifstream infile_mesh(vtkfile.c_str());

            // First, read the mesh data
            mesh.read_vtk(infile_mesh);

            std::ifstream infile(vtkfile.c_str());

            // Read point data
            std::istringstream ss;
            std::string str, name, type, line;
            int numcomp;
            unsigned long int Npoints;

            while(!infile.eof())
            {
                std::getline( infile, line);

                if( line[0] == 'P')
                {
                    // Read the header and determine the number of points
                    if( line.size() > 9 && line.substr(0,10) == "POINT_DATA")
                    {
                        ss.clear();
                        ss.str(line);
                        ss >> str >> Npoints;
                    }
                }
                else if( line[0] == 'S')
                {
                    // Read the scalar field data
                    if( line.size() > 6 && line.substr(0,7) == "SCALARS")
                    {
                        ss.clear();
                        ss.str(line);
                        ss >> str >> name >> type >> numcomp;

                        // Read LOOKUP_TABLE line
                        std::getline( infile, line);

                        // Read data
                        std::cout << "begin reading data" << std::endl;
                        std::vector<double> data(Npoints);
                        for( unsigned int i=0; i<Npoints; i++)
                        {
                            infile >> data[i];
                        }
                        std::cout << "  done" << std::endl;

                        // Construct the scalar field as a PField object:
                        //     Pass the field name, mesh, data, and field type
                        std::vector<std::string> var_name(DIM);
                        std::vector<std::string> var_description(DIM);

                        if( DIM >= 2)
                        {
                            var_name[0] = "x";
                            var_description[0] = "x coordinate";
                            var_name[1] = "y";
                            var_description[1] = "y coordinate";
                        }
                        if( DIM >= 3)
                        {
                            var_name[2] = "z";
                            var_description[2] = "z coordinate";
                        }

                        std::cout << "Construct PField '" << name << "'" << std::endl;
                        scalar_field.push_back( PField<Coordinate, double, DIM>( name, var_name, var_description, mesh, data, 0.0) );
                        std::cout << "  done" << std::endl;
                    }
                }
                // Alternative field descriptor used by ParaView (holds the same 
                // information as the "SCALAR" line above)
                else if( line[0] == 'F')
                {
                    if( line.size() > 14 && line.substr(0,15) == "FIELD FieldData")
                    {
                        ss.clear();
                        ss.str(line);
                        ss >> str >> numcomp;

                        // Read LOOKUP_TABLE line
                        std::getline( infile, line);

                        ss.clear();
                        ss.str(line);
                        ss >> name >> numcomp >> Npoints >> type;

                        // Read the field data
                        std::cout << "begin reading data" << std::endl;
                        std::vector<double> data(Npoints);
                        for( unsigned int i=0; i<Npoints; i++)
                        {
                            infile >> data[i];
                        }
                        std::cout << "  done" << std::endl;

                        // Construct the scalar field as a PField object:
                        //     Pass the field name, mesh, data, and field type
                        std::vector<std::string> var_name(DIM);
                        std::vector<std::string> var_description(DIM);

                        if( DIM >= 2)
                        {
                            var_name[0] = "x";
                            var_description[0] = "x coordinate";
                            var_name[1] = "y";
                            var_description[1] = "y coordinate";
                        }
                        if( DIM >= 3)
                        {
                            var_name[2] = "z";
                            var_description[2] = "z coordinate";
                        }

                        std::cout << "Construct PField '" << name << "'" << std::endl;
                        scalar_field.push_back( PField<Coordinate, double, DIM>( name, var_name, var_description, mesh, data, 0.0) );
                        std::cout << "  done" << std::endl;

                    }
                }
            }
            infile.close();
        }

        // Function for returning scalar field data
        PField<Coordinate, double, DIM>& find_scalar_field(std::string name)
        {
            for( unsigned int i=0; i<scalar_field.size(); i++)
            {
                if( scalar_field[i].name() == name)
                    return scalar_field[i];
            }
            throw std::invalid_argument("Could not find scalar_field named '" + name + "'" );
        }
    };
}

#endif
