void read_baron_results(double* objval, double* sol, int sol_size)
 {

     //first read .time file
     std::ifstream timefile("tim.lst");
     if (!timefile)
     {
         std::cerr << "Error: could not open file for writing.\n";
         Assert(1 == 0);
     }
     std::ostringstream buffer;
     buffer << timefile.rdbuf(); // read full file content into buffer
     std::string content = buffer.str();

     // Remove trailing newline(s)
     content.erase(std::find_if(content.rbegin(), content.rend(),
                                [](unsigned char ch) { return !std::isspace(ch); }).base(), content.end());

     std::istringstream iss(content);
     std::vector<std::string> tokens;
     std::string s;
     while (iss >> s) {
         tokens.push_back(s);
     }

     //We always have a minimization problem
     double dual_objval = std::stod(tokens[5]);
     *objval = std::stod(tokens[6]);

     //for max problem
     /*
       double dual_objval = std::stod(tokens[6]);
       *objval = std::stod(tokens[5]);

      */

     int solver_status = std::stoi(tokens[7]);
     int moedl_status = std::stoi(tokens[8]);
     int nodeopt = std::stoi(tokens[11]);

     if(nodeopt != -3)
     {
         std::ifstream resultfile("res.lst");
         if (!resultfile)
         {
             std::cerr << "Error: could not open file for writing.\n";
             Assert(1 == 0);
         }

         std::string line;
         bool32 found = false;

         // Loop through lines to find the target
         while (std::getline(resultfile, line))
         {
             if (line.rfind("The best solution found", 0) == 0) {  // startswith check
                 found = true;
                 break;
             }
         }

         if (!found)
         {
             std::cerr << "Error: could not find solution in result file.\n";
             Assert(1 == 0);
         }

         std::getline(resultfile, line); // discard
         std::getline(resultfile, line);

         std::regex number_regex(R"(\d+)");
         while (std::getline(resultfile, line)) {
             // Trim right-side whitespace (similar to chomp in Julia)
             line.erase(line.find_last_not_of(" \t\r\n") + 1);

             // If line is empty, break the loop
             if (line.empty()) break;

             // Split line by whitespace
             std::istringstream result_iss(line);
             std::vector<std::string> parts;
             std::string part;
             while (result_iss >> part) {
                 parts.push_back(part);
             }

             // Check if we have at least 4 fields
             if (parts.size() < 4) {
                 std::cerr << "Expected at least 4 parts per line, got: " << line;
                 Assert(1 == 0);
             }

             // Match digits in parts[0]
             std::smatch match;
             if (!std::regex_search(parts[0], match, number_regex)) {
                 std::cerr << "Cannot find appropriate variable index from " << parts[0];
                 Assert(1 == 0);
             }

             int v_idx = std::stoi(match[0]);
             double v_val = std::stod(parts[2]);

             sol[v_idx] = v_val;
         }

     }
     
 }








halfspace minimum_halfspace_slices(Array<Matrix<Rational>> slices_vertices,  Matrix<Rational> all_vertices,  Array<Array<Set<Int>>> slices_triangulations, Array<Array<Set<Int>>> slices_edges,  Rational TotalVolume,  Array<Rational> slices_volumes,  Array<Int> slices_dimensions,int64 ambient_dimension, int centerpoint_slice_index, Vector<Rational> centerpoint, int iteration)
{


    int num_slices = slices_vertices.size();
   
    int gurobi_infeasibilities = 0;

    Matrix<Rational> normals(all_vertices.rows(), all_vertices.cols());  // Create an empty matrix with correct size
    for (int i = 0; i < all_vertices.rows(); i++)
    {
        normals[i] = all_vertices.row(i) - centerpoint;  // Assign the same vector to each row
    }
   

    BigObject HA("HyperplaneArrangement<Rational>");

    HA.take("HYPERPLANES") << normals;

    IncidenceMatrix<> signatures = HA.give("CHAMBER_SIGNATURES");

    BigObject fan("PolyhedralFan<Rational>");
    fan = HA.give("CHAMBER_DECOMPOSITION");
    Matrix<Rational> all_rays = fan.give("RAYS");

    
    //NOTE: Indices refer to RAYS.
    IncidenceMatrix<> maximal_cones = fan.give("MAXIMAL_CONES");
      
    //Loop over all signaurtes
    cout << signatures << endl;

   
   
    std::vector<int64> minimal_volume_signature_indices;
    std::vector<real64> minimal_volumes;
    int row_dimension = all_vertices.cols();
   
    Matrix<Rational> minimal_directions = zero_matrix<Rational>(signatures.rows(), row_dimension);
    
    
    for(int i = 0; i < signatures.rows(); i++)
    {
        cout << "-------------------------------------------------------------------------------------" << endl;
        cout << "Computing Halfspace Depth for Vertices: ";

        Matrix<Rational> rays;
        //Vector<Rational> defining_normal = zero_vector<Rational>(all_rays.cols());
        //walks through all rays
        for(int j = 0; j < maximal_cones.cols(); j++)
        {
            if(maximal_cones(i,j))
            {
                rays /= all_rays.row(j);
                // defining_normal += all_rays.row(j);
            }
        }

        
        //NOTE: We can't use the rays as inequality here. We need to convert it to a h-representation!
        BigObject cone("Cone<Rational>");
        cone.take("INPUT_RAYS") << rays;
        Matrix<Rational> cone_h_representation = cone.give("FACETS");
        Matrix<Rational> cone_h_lin = cone.give("LINEAR_SPAN");

        
     

       
        //NOTE: (Update): This is left out for now:
        // The problem was that we allowed a direction of 1 in one direction
        // So with lower dimensional one slices we alwas get a unit vector 
        // We later killed that direction, so no progress was possible

               //TODO:Generalize for arbitrary dim
               //TODO: Should I deal with this case ?
               // if num Slices is 1 the slice should be profiede full dimensional
               //NOTE: THIS IS NEED: Otherwise infinite loop, in feasability checking
               // of gradient descent

               /*
        bool32 all_in_one_slice = true;
        for(int j = 0; j < signatures.cols(); j++)
        {
            if(signatures(i,j))
            {
                if(all_vertices[j][0] != centerpoint[0])
                {
                    all_in_one_slice = false;
                }
            }
        }
           
        if(all_in_one_slice && (num_slices == 1))
        {
            cout << "rays " << endl;
            cout << rays << endl;
            
            cout << "cone_h_representation" << endl;
            cout << cone_h_representation << endl;

            
            Vector<Rational> linear_span_cone_vec = {1, 0, 0};
            cone_h_representation /= linear_span_cone_vec;
            cone_h_representation /= -linear_span_cone_vec;
            

            //Vector<Rational> linear_span_cone_vec = {1, 100, 0};
            //cone_h_representation /= linear_span_cone_vec;
            
        }
               */
           
        

        cout << "rays " << endl;
        cout << rays << endl;
            
        cout << "cone_h_representation" << endl;
        cout << convert_to<double>(cone_h_representation) << endl;


        //RESCALE cone_h_representation
        //TODO: BETTER FOR BASE TWO exponents?
        real64 highest_exponent = 1.0;
        
        for(int j = 0; j < cone_h_representation.rows(); j++)
        {
            for(int k = 0; k < cone_h_representation.cols(); k++)
            {
                real64 exponent = std::floor(std::log10(std::fabs(convert_to<double>(cone_h_representation(j,k)))));
                if(exponent > highest_exponent)
                {
                    highest_exponent = exponent;
                }
                    
            }
        }

        int scaling_exponent = static_cast<int>(highest_exponent)-1;
        
        
        double tentimes = std::pow(10.0, scaling_exponent);
        
        if(tentimes > 0)
        {
            cone_h_representation /= convert_to<Rational>(tentimes);
        }

        cout << "cone_h_representation (scaled)" << endl;
        cout << convert_to<double>(cone_h_representation) << endl;


        //THIS IS BETTER

        BigObject newchamber("Cone<Rational>");
        Matrix<Rational> inequalities(normals);
        inequalities.minor(~signatures[i], All) *= -1;
        newchamber.take("INEQUALITIES") << inequalities;
        Matrix<Rational> cone_h_rep = newchamber.give("FACETS");



        cout << "new_cone_h_rep" << endl;
        cout << convert_to<double>(cone_h_rep) << endl;

        cone_h_representation = cone_h_rep;

        

    
        
        Matrix<Rational> sub_vertices;
        Matrix<Rational> conter_vertices;
        
       
        //TODO: Is there a more efficient methode for this than cunstructing a set?
        Set<Int> sub_points;
        Set<Int> conter_points;
        std::unordered_map<int, int>  sub_points_array = {};
        std::unordered_map<int, int>  conter_points_array = {};
        int64 sub_index = 0;
        int64 conter_index = 0;
        for(int j = 0; j < signatures.cols(); j++)
        {
            
            if(signatures(i,j))
            {
                cout << j << ", ";
                sub_vertices /= all_vertices.row(j);
                sub_points.collect(j);
                sub_points_array.insert({sub_index, j});
                sub_index++;
            }
            else
            {
                conter_vertices /= all_vertices.row(j);
                conter_points.collect(j);
                conter_points_array.insert({conter_index, j});
                conter_index++;
            }
        }



       



        /*_________sub_polytope___________*/


        

        /*_________conter_polytope___________*/
        

        
        int32 number_lambda = 0;

        
        Set<Int> sub_boundary_points;
        Set<Int> conter_boundary_points;

        Matrix<Rational> sub_boundary_vertices;
        Matrix<Rational> conter_boundary_vertices;

        std::vector<long> sub_boundary_array = {};
        std::vector<long> conter_boundary_array = {};

        Matrix<Int> edge_map = -1*ones_matrix<Int>(all_vertices.rows(), all_vertices.rows());

        
        int vertices_so_far = 0;
     
        for(int j = 0; j < num_slices; j++)
        {

            Array<Set<Int>> slice_edges = slices_edges[j];
            for(int k = 0; k < slice_edges.size(); k++)
            {
                if((signatures(i,vertices_so_far + slice_edges[k].front())
                    && !signatures(i,vertices_so_far + slice_edges[k].back())))
                {
                
                    sub_boundary_vertices /= all_vertices.row(vertices_so_far + slice_edges[k].front());
                    sub_boundary_points.collect(vertices_so_far + slice_edges[k].front());
                    sub_boundary_array.push_back(vertices_so_far + slice_edges[k].front());

                    conter_boundary_vertices /= all_vertices.row(vertices_so_far + slice_edges[k].back());
                    conter_boundary_points.collect(vertices_so_far + slice_edges[k].back());
                    conter_boundary_array.push_back(vertices_so_far + slice_edges[k].back());

                    edge_map(vertices_so_far + slice_edges[k].front(),
                             vertices_so_far + slice_edges[k].back()) = number_lambda;
                    number_lambda++;

                
                
                }
                else if((signatures(i,vertices_so_far + slice_edges[k].back())
                         && !signatures(i,vertices_so_far + slice_edges[k].front())))
                {
                    sub_boundary_vertices /= all_vertices.row(vertices_so_far + slice_edges[k].back());
                    sub_boundary_points.collect(vertices_so_far + slice_edges[k].back());
                    sub_boundary_array.push_back(vertices_so_far + slice_edges[k].back());
              
                
                    conter_boundary_vertices /= all_vertices.row(vertices_so_far + slice_edges[k].front());
                    conter_boundary_points.collect(vertices_so_far + slice_edges[k].front());
                    conter_boundary_array.push_back(vertices_so_far + slice_edges[k].front());

                    edge_map(vertices_so_far + slice_edges[k].back(),
                             vertices_so_far + slice_edges[k].front()) = number_lambda;
                    number_lambda++;
                
                }
            }

             vertices_so_far += slices_vertices[j].rows();
        }

        
       
        
        Vector<Vector<Matrix<Rational>>> quartobj_coeff(number_lambda);
        for(int w = 0; w < number_lambda; w++)
        {
            quartobj_coeff[w] = Vector<Matrix<Rational>>(number_lambda);
            for(int v = 0; v < number_lambda; v++)
            {
                quartobj_coeff[w][v] = zero_matrix<Rational>(number_lambda, number_lambda);
            }
        }
        Vector<Matrix<Rational>> triobj_coeff(number_lambda);
        for(int w = 0; w < number_lambda; w++)
        {
            triobj_coeff[w] = zero_matrix<Rational>(number_lambda, number_lambda);
        }
        Matrix<Rational> qobj_coeff = zero_matrix<Rational>(number_lambda, number_lambda);
        Vector<Rational> lobj_coeff = zero_vector<Rational>(number_lambda);
        Rational FixVolume = 0;
        Rational ConterVolume = 0;
        vertices_so_far = 0;
        
        for(int j = 0; j < num_slices; j++)
        {
            Array<Set<Int>> slice_triangulation = slices_triangulations[j];

            


            for(int k = 0; k < slice_triangulation.size(); k++)
            {

                std::vector<lambda_u_constraint> constraints;
                Vector<Int> triangle = Vector<Int>(slice_triangulation[k]);
                //compute the total volume of the triangle:

                cout <<"slice " << j << " triag: " <<  triangle << endl;

                for(int l = 0; l < triangle.size(); l++)
                {
                    if(sub_boundary_points.contains(vertices_so_far + triangle[l]))
                    {
                        for(int m = 0; m < triangle.size(); m++)
                        {
                            if(conter_boundary_points.contains(vertices_so_far + triangle[m]))
                            {
                                lambda_u_constraint cutting_edge = {};
                                cutting_edge.start_vec_index = vertices_so_far + triangle[l];
                                cutting_edge.end_vec_index = vertices_so_far + triangle[m];
                                constraints.push_back(cutting_edge);
                           
                            }

                        }
                       
                    }
                
                }

                
                //NOTE: HERE THE TRIANGULATION IS A BOUNDARY TRIANGULATION

                switch (slices_dimensions[j])
                {
                    case 1:
                        // We Always have 0 or 1 contraits
                        if(j == centerpoint_slice_index)
                        {
                            if(constraints.size() > 0)
                            {
                                Assert(constraints.size() == 1);
                                Vector<Rational> vec0 = all_vertices.row(constraints[0].start_vec_index);
                                Vector<Rational> vec1 = all_vertices.row(constraints[0].end_vec_index);

                                

                                //NOTE: THIS IS a slice form polymake!!
                                vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                Vector<Rational> proj_centerpoint = centerpoint.slice(range_from(ambient_dimension - slices_dimensions[j]));

                                Assert(vec0.dim() == 1);

                                FixVolume += AbsoluteValue(proj_centerpoint[0] - vec0[0]);
                                ConterVolume += AbsoluteValue(vec1[0] - proj_centerpoint[0]);
 
                            }
                            else
                            {
                                //This case should never happen!!
                                //Maybe if the centerpoint is the vertex...
                                Assert(1 == 0);
                            }
                        }
                        else
                        {
                            if(constraints.size() > 0)
                            {
                                Assert(constraints.size() == 1);
                                Vector<Rational> vec0 = all_vertices.row(constraints[0].start_vec_index);
                                Vector<Rational> vec1 = all_vertices.row(constraints[0].end_vec_index);

                                //NOTE: THIS IS a slice form polymake!!
                                vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));

                                Assert(vec0.dim() == 1);

                                Rational lambda1 = vec1[0] - vec0[0];

                                lobj_coeff[edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index)]
                                    += AbsoluteValue(lambda1);
                            }
                            else
                            {
                                Vector<Rational> vec0 = all_vertices.row(vertices_so_far + triangle[0]);
                                Vector<Rational> vec1 = all_vertices.row(vertices_so_far + triangle[1]);


                                //NOTE: THIS IS a slice form polymake!!
                                vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));

                                Assert(vec0.dim() == 1);

                                Rational Volume = vec1[0] - vec0[0];

                                if(sub_points.contains(vertices_so_far + triangle[0]))
                                {
                                    FixVolume += AbsoluteValue(Volume);
                                }
                                else
                                {
                                    ConterVolume += AbsoluteValue(Volume);
                                }
                            }
                            
                        }
                        break;
                            
                    case 2:
                        // We Always have 0 or 1 contraits

                        if(j == centerpoint_slice_index)
                        {
                            if(constraints.size() > 0)
                            {
                                Assert(constraints.size() == 1);

                                Vector<Rational> vec0 = all_vertices.row(constraints[0].start_vec_index);
                                Vector<Rational> vec1 = all_vertices.row(constraints[0].end_vec_index);

                                //NOTE: THIS IS a slice form polymake!!
                                vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                Vector<Rational> proj_centerpoint = centerpoint.slice(range_from(ambient_dimension - slices_dimensions[j]));

                                Rational lambda_1 = det(vector2row(vec1-vec0) / (proj_centerpoint - vec0));

                                lobj_coeff[edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index)]
                                    += Rational(1,2)*AbsoluteValue(lambda_1);

                            }
                            else
                            {

                                Vector<Rational> vec0 = all_vertices.row(vertices_so_far + triangle[0]);
                                Vector<Rational> vec1 = all_vertices.row(vertices_so_far + triangle[1]);


                                //NOTE: THIS IS a slice form polymake!!
                                vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                Vector<Rational> proj_centerpoint = centerpoint.slice(range_from(ambient_dimension - slices_dimensions[j]));

                                Rational Volume =  det(vector2row(vec1-vec0) / (proj_centerpoint - vec0));

                                if(sub_points.contains(vertices_so_far + triangle[0]))
                                {
                                    FixVolume += Rational(1,2)*AbsoluteValue(Volume);
                                }
                                else
                                {
                                    ConterVolume += Rational(1,2)*AbsoluteValue(Volume);
                                }
                                
                            }
                        }
                        else
                        {
                            if(constraints.size() > 0)
                            {
                                Assert(constraints.size() == 2);

                                if(constraints[0].start_vec_index == constraints[1].start_vec_index)
                                {
                                    /* case:
                                       current |----------| current_conter
                                       | \        |
                                       |   \      |
                                       |      \   |
                                       |        \ |
                                       next |----------| next_conter

                                    */

                        


                                    Vector<Rational> vec0 = all_vertices.row(constraints[0].start_vec_index);
                                    Vector<Rational> vec1 = all_vertices.row(constraints[0].end_vec_index);
                                    Vector<Rational> vec2 = all_vertices.row(constraints[1].end_vec_index);


                                    //NOTE: THIS IS a slice form polymake!!
                                    vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));


                                    //Computing Volume via determinant
                                    Rational lambda_1_lambda_2 = det(vector2row(vec1-vec0) / (vec2-vec0));

                                    cout << "number_lambda "<< number_lambda << endl;
                                    cout << "l1 " << constraints[0].start_vec_index << endl;
                                    cout << "l2 " << constraints[0].end_vec_index << endl;
                                    cout << "l3 " << constraints[1].end_vec_index << endl;
                        
                                    qobj_coeff(edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index),
                                               edge_map(constraints[1].start_vec_index, constraints[1].end_vec_index)) =
                                        Rational(1,2)*AbsoluteValue(lambda_1_lambda_2);
                               
                                

                                }

                                else
                                {

                                    /* case:
                                       current |----------| current_conter
                                       |        / |
                                       |      /   |
                                       |   /      |
                                       | /        |
                                       next |----------| next_conter

                                    */

                      

                                    //Compute Volume:
                    
                                    //Triangle 1: 
                                        Vector<Rational> vec0 = all_vertices.row(constraints[0].start_vec_index);
                                    Vector<Rational> vec1 = all_vertices.row(constraints[1].start_vec_index);
                                    Vector<Rational> vec2 = all_vertices.row(constraints[0].end_vec_index);


                                    //NOTE: THIS IS a slice form polymake!!
                                    vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));



                                    //Computing Volume via determinant
                                    Rational lambda_1_lambda_2 = det(vector2row(vec2-vec0) / (vec2-vec1));
                                    Rational lambda_1 = det(vector2row(vec2-vec0) / (vec1-vec0));
                                            

                                    //Triangle 2: 
                                    //Computing Volume via determinant
                       
                                    Rational lambda_2 = det(vector2row(vec2-vec1) / (vec1-vec0));
                                    Rational constant = det(vector2row(vec1-vec0) / (vec1-vec0));
                                        
                                    Assert(convert_to<double>(constant) == 0);

                                    if(det(vector2row(vec2-vec0) / (vec1-vec0)) >= 0)
                                    {
                                        //NOTE in lobj_coeff += necessary since we might push linear terms ahead

                                        qobj_coeff(edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index),
                                                   edge_map(constraints[1].start_vec_index, constraints[1].end_vec_index)) =
                                            Rational(1,2)*lambda_1_lambda_2;
                             
                                        lobj_coeff[edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index)]
                                            += Rational(1,2)*lambda_1;

                           
                                        lobj_coeff[edge_map(constraints[1].start_vec_index, constraints[1].end_vec_index)]
                                            += Rational(1,2)*lambda_2;

                                    }

                                    else
                                    {
                                        qobj_coeff(edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index),
                                                   edge_map(constraints[1].start_vec_index, constraints[1].end_vec_index)) =
                                            Rational(1,2)*(-lambda_1_lambda_2);
                             
                                        lobj_coeff[edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index)]
                                            += Rational(1,2)*(-lambda_1);

                           
                                        lobj_coeff[edge_map(constraints[1].start_vec_index, constraints[1].end_vec_index)]
                                            += Rational(1,2)*(-lambda_2);
                                    }

                        
                            
                                }
                            }
                            else
                            {

                                Vector<Rational> vec0 = all_vertices.row(vertices_so_far + triangle[0]);
                                Vector<Rational> vec1 = all_vertices.row(vertices_so_far + triangle[1]);
                                Vector<Rational> vec2 = all_vertices.row(vertices_so_far + triangle[2]);


                                //NOTE: THIS IS a slice form polymake!!
                                vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                   
                                Rational Volume =  det(vector2row(vec1-vec0) / (vec2 - vec0));

                                if(sub_points.contains(vertices_so_far + triangle[0]))
                                {
                                    FixVolume += Rational(1,2)*AbsoluteValue(Volume);
                                }
                                else
                                {
                                    ConterVolume += Rational(1,2)*AbsoluteValue(Volume);
                                }
                                
                            }

                        }
                        break;
                    case 3:

                        if(j == centerpoint_slice_index)
                        {
                            // We Always have 0 or 2 contraits
                            if(constraints.size() > 0)
                            {
                                Assert(constraints.size() == 2);
                                if(constraints[0].start_vec_index == constraints[1].start_vec_index)
                                {
                                    /* case:
                                       current |----------| current_conter
                                       | \        |
                                       |   \      |
                                       |      \   |
                                       |        \ |
                                       next |----------| next_conter

                                    */

                        


                                    Vector<Rational> vec0 = all_vertices.row(constraints[0].start_vec_index);
                                    Vector<Rational> vec1 = all_vertices.row(constraints[0].end_vec_index);
                                    Vector<Rational> vec2 = all_vertices.row(constraints[1].end_vec_index);

                                    vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    Vector<Rational> proj_centerpoint = centerpoint.slice(range_from(ambient_dimension - slices_dimensions[j]));


                                    //Computing Volume via determinant
                                    Rational lambda_1_lambda_2 = cross_product3(vec1-vec0, vec2-vec0)*(proj_centerpoint-vec0);

                        
                                    qobj_coeff(edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index),
                                               edge_map(constraints[1].start_vec_index, constraints[1].end_vec_index)) =
                                        Rational(1,6)*AbsoluteValue(lambda_1_lambda_2);
                               
                                
                                    Vector<Rational> l1_nominator= proj_centerpoint-vec0;
                                    Vector<Rational> l1_denominator = vec1-vec0;

                                }

                                else
                                {

                                    /* case:
                                       current |----------| current_conter
                                       |        / |
                                       |      /   |
                                       |   /      |
                                       | /        |
                                       next |----------| next_conter

                                    */

                      

                                    //Compute Volume:
                    
                                    //Triangle 1: 
                                        Vector<Rational> vec0 = all_vertices.row(constraints[0].start_vec_index);
                                    Vector<Rational> vec1 = all_vertices.row(constraints[1].start_vec_index);
                                    Vector<Rational> vec2 = all_vertices.row(constraints[0].end_vec_index);

                                    vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    Vector<Rational> proj_centerpoint = centerpoint.slice(range_from(ambient_dimension - slices_dimensions[j]));


                                    //Computing Volume via determinant
                                    Rational lambda_1_lambda_2 = cross_product3(vec2-vec0, vec2-vec1)*(proj_centerpoint-vec0);
                                    Rational lambda_1 =cross_product3(vec2-vec0, vec1-vec0)*(proj_centerpoint-vec0);

                                    //Triangle 2: 
                                    //Computing Volume via determinant
                       
                                    Rational lambda_2 = cross_product3(vec2-vec1, vec1-vec0)*(proj_centerpoint-vec0);
                                    Rational constant = cross_product3(vec1-vec0, vec1-vec0)*(proj_centerpoint-vec0);
                                    Assert(convert_to<double>(constant) == 0);

                                    if((cross_product3(vec2-vec0, vec1-vec0)*(proj_centerpoint - vec0)) >= 0)
                                    {
                                        //NOTE in lobj_coeff += necessary since we might push linear terms ahead

                                        qobj_coeff(edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index),
                                                   edge_map(constraints[1].start_vec_index, constraints[1].end_vec_index)) =
                                            Rational(1,6)*lambda_1_lambda_2;
                             
                                        lobj_coeff[edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index)]
                                            += Rational(1,6)*lambda_1;

                           
                                        lobj_coeff[edge_map(constraints[1].start_vec_index, constraints[1].end_vec_index)]
                                            += Rational(1,6)*lambda_2;

                                    }

                                    else
                                    {
                                        qobj_coeff(edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index),
                                                   edge_map(constraints[1].start_vec_index, constraints[1].end_vec_index)) =
                                            Rational(1,6)*(-lambda_1_lambda_2);
                             
                                        lobj_coeff[edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index)]
                                            += Rational(1,6)*(-lambda_1);

                           
                                        lobj_coeff[edge_map(constraints[1].start_vec_index, constraints[1].end_vec_index)]
                                            += Rational(1,6)*(-lambda_2);
                                    }

                        
                            
                                }
                    
                            }

                            else
                            {
                                Vector<Rational> vec0 = all_vertices.row(vertices_so_far + triangle[0]);
                                Vector<Rational> vec1 = all_vertices.row(vertices_so_far + triangle[1]);
                                Vector<Rational> vec2 = all_vertices.row(vertices_so_far + triangle[2]);

                                vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                Vector<Rational> proj_centerpoint = centerpoint.slice(range_from(ambient_dimension - slices_dimensions[j]));

                                Rational Volume = cross_product3(vec1-vec0, vec2-vec0)*(proj_centerpoint-vec0);

                                if(sub_points.contains(vertices_so_far + triangle[0]))
                                {
                                    FixVolume += Rational(1,6)*AbsoluteValue(Volume);
                                }
                                else
                                {
                                    ConterVolume += Rational(1,6)*AbsoluteValue(Volume);
                                }

                            }
                        }
                        else
                        {
                            // We Always have 0, 3 or 4 contraits
                            if(constraints.size() == 4)
                            {
                                //split vertices 2 by 2

                                /*

                                  1  ----   4
                                  | 2  --   | 5
                                  |/        |/
                                  0   ----  3

                                  NOTE: vec1 == vec4 and vec2 == vec5
                                */

                                int64 second_constraint_index = 0;

                                
                                for(int n = 1; n < constraints.size(); n++)
                                {
                                    if(constraints[0].start_vec_index != constraints[n].start_vec_index)
                                    {
                                        second_constraint_index = n;
                                    }

                                }


                                Vector<Rational> vec0 = all_vertices.row(constraints[0].start_vec_index);
                                Vector<Rational> vec1 = all_vertices.row(constraints[0].end_vec_index);
                               
                                int64 end_vec_index;
                                   
                                int64 add_edges = 0;
                               
                                for(int n = 1; n < constraints.size(); n++)
                                {
                                    if(constraints[0].start_vec_index == constraints[n].start_vec_index)
                                    {
                                        add_edges++;
                                        end_vec_index = n;
                                        break;
                                    }
                                }
                                Assert(add_edges == 1);

                                Vector<Rational> vec2 = all_vertices.row(end_vec_index);
                            
                                        
                                Vector<Rational> vec3 = all_vertices.row(constraints[second_constraint_index].start_vec_index);
                                Vector<Rational> vec4 = all_vertices.row(constraints[0].end_vec_index);
                                Vector<Rational> vec5 = all_vertices.row(end_vec_index);


                                //NOTE: THIS IS a slice form polymake!!
                                vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec3 = vec3.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec4 = vec4.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec5 = vec5.slice(range_from(ambient_dimension - slices_dimensions[j]));

                                int64 l1 = edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index);
                                int64 l2 = edge_map(constraints[0].start_vec_index, end_vec_index);
                                int64 l4 = edge_map(constraints[second_constraint_index].start_vec_index, constraints[0].end_vec_index);

                                int64 l5 = edge_map(constraints[second_constraint_index].start_vec_index, end_vec_index);
                                        
                                        
                                //Triangulation
                                //{0,1,2,3},{1,2,3,4},{2,3,4,5}

                                //{0,1,2,3}
                                //determinant with base point 0 ---> orientation "-"
                                //det(l1*(vec1 - vec0), l2*(vec2-vec0), vec3 - vec0)
                                Rational lambda_1_lambda_2 =  det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec3 - vec0));
                               

                                //{1,2,3,4}
                                //determinant with base point 3 ---> orientation "-"
                                //det(l1*(vec1 - vec0) - (vec3 - vec0), l2*(vec2-vec0)-
                                //    (vec3- vec0), l4(vec4 - vec 3))
                                        
                                       
                                Rational  lambda_1_lambda_2_lambda_4 =  det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec4 - vec3));

                                //0: det(vec1 -vec 0, vec3-vec0, vec1-vec0) all in one plane!
                                Rational  lambda_1_lambda_4 =  det(vector2row(vec1-vec0) / -(vec3 - vec0) / (vec4 - vec3));
                                        
                                Rational  lambda_2_lambda_4 =  det(vector2row(-(vec3-vec0)) / (vec2 - vec0) / (vec4 - vec3));

                                //0:
                                Rational  lambda_4 =  det(vector2row(-(vec3-vec0)) / -(vec3 - vec0) / (vec3 - vec0));

                                        


                                //{2,3,4,5}
                                //determiant with basepoint 3 --> orientation "+"
                                /*
                                  det(l2*(vec2 - vec0) - (vec3-vec0), l4*(vec4-vec3),
                                  l5*(vec5 - vec3))
                                */
                                Rational  lambda_2_lambda_4_lambda_5 =  det(vector2row(vec2-vec0) / (vec4 - vec3) / (vec5 - vec3));


                                Rational  lambda_4_lambda_5 =  det(vector2row(-(vec3-vec0)) / (vec4 - vec3) / (vec5 - vec3));

                                
                                //NOTE: Determin sign based on first triangle, this on is non singular since it arises from the triagle!
                                if(det(vector2row(vec1-vec0) / (vec2-vec0) / (vec3-vec0)) > 0)  
                                {
                                    //1.
                                    qobj_coeff(l1,l2) += Rational(1,6)*lambda_1_lambda_2;
                                    

                                    //2.
                                    triobj_coeff[l1](l2, l4) +=  Rational(1,6)*lambda_1_lambda_2_lambda_4;
                                    qobj_coeff(l1,l4) += Rational(1,6)*lambda_1_lambda_4;

                                    qobj_coeff(l2,l4) += Rational(1,6)*lambda_2_lambda_4;

                                    lobj_coeff[l4] +=  Rational(1,6)*lambda_4;


                                    //3.
                                    triobj_coeff[l2](l4, l5) +=  Rational(1,6)*-lambda_2_lambda_4_lambda_5;
                                    qobj_coeff(l4,l5) += Rational(1,6)*-lambda_4_lambda_5;
                                                       
                                }
                                else
                                {
                                     //1.
                                    qobj_coeff(l1,l2) += Rational(1,6)*-lambda_1_lambda_2;
                                    

                                    //2.
                                    triobj_coeff[l1](l2, l4) +=  Rational(1,6)*-lambda_1_lambda_2_lambda_4;
                                    qobj_coeff(l1,l4) += Rational(1,6)*-lambda_1_lambda_4;

                                    qobj_coeff(l2,l4) += Rational(1,6)*-lambda_2_lambda_4;

                                    lobj_coeff[l4] +=  Rational(1,6)*-lambda_4;


                                    //3.
                                    triobj_coeff[l2](l4, l5) +=  Rational(1,6)*lambda_2_lambda_4_lambda_5;
                                    qobj_coeff(l4,l5) += Rational(1,6)*lambda_4_lambda_5;

                                }

                               

                                
                            }
                            else if(constraints.size() == 3)
                            {
                                if(constraints[0].start_vec_index == constraints[1].start_vec_index)
                                {
                                    /*

                                      1  
                                      | 2 
                                      |/        
                                      0--3

                                    */

                                    Vector<Rational> vec0 = all_vertices.row(constraints[0].start_vec_index);
                                    Vector<Rational> vec1 = all_vertices.row(constraints[0].end_vec_index);
                                    Vector<Rational> vec2 = all_vertices.row(constraints[1].end_vec_index);
                                    Vector<Rational> vec3 = all_vertices.row(constraints[2].end_vec_index);
                                    vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec3 = vec3.slice(range_from(ambient_dimension - slices_dimensions[j]));

                                    int64 l1 = edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index);
                                    int64 l2 = edge_map(constraints[1].start_vec_index, constraints[1].end_vec_index);
                                    int64 l3 = edge_map(constraints[2].start_vec_index, constraints[2].end_vec_index);
                                    Rational lambda_1_lambda_2_lambda_3 =  det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec3 - vec0));

                                    triobj_coeff[l1](l2, l3) +=  Rational(1,6)*AbsoluteValue(lambda_1_lambda_2_lambda_3);

                                }

                                else
                                {
                                     /*
                                       2 -
                                      | \  \
                                      |  |  \
                                      |  |   -
                                      |   1   \
                                      5   / \   -
                                      | 4    \   \
                                      |/       \  \
                                      x-----3-----0

                                    */

                                     //Triangulation
                                    //{0,1,2,3},{1,2,3,4},{2,3,4,5}
                                    //NOTE: WE have the same triangulation, but diffrent cutting edges

                                    Vector<Rational> vec0 = all_vertices.row(constraints[0].start_vec_index);
                                    Vector<Rational> vec1 = all_vertices.row(constraints[1].start_vec_index);
                                    Vector<Rational> vec2 = all_vertices.row(constraints[2].start_vec_index);
                                    Vector<Rational> vec3 = all_vertices.row(constraints[0].end_vec_index);
                                    Vector<Rational> vec4 = all_vertices.row(constraints[1].end_vec_index);
                                    Vector<Rational> vec5 = all_vertices.row(constraints[2].end_vec_index);

                                     //NOTE: vec3 == vec4 == vec5

                                    vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec3 = vec3.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec4 = vec4.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec5 = vec5.slice(range_from(ambient_dimension - slices_dimensions[j]));

                                    int64 l1 = edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index);
                                    int64 l2 = edge_map(constraints[1].start_vec_index, constraints[1].end_vec_index);
                                    int64 l3 = edge_map(constraints[2].start_vec_index, constraints[2].end_vec_index);

                                    
                                    //{0,1,2,3}
                                    //determinant with base point 0 --> orientation "-"
                                    //det((vec1 - vec0), (vec2-vec0), l1*(vec3 - vec0))
                                    Rational lambda_1=  det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec3 - vec0));

                                    //{1,2,3,4}
                                    //determinant with base point 1 --> orientation "-"
                                    //det((vec2 - vec1), l1*(vec3 - vec0) + (vec0 -vec1), l2*(vec4-vec1))
                                        
                                       
                                    Rational  lambda_1_lambda_2 =  det(vector2row(vec2-vec1) / (vec3 - vec0) / (vec4 - vec1));

                                    Rational  lambda_2 =   det(vector2row(vec2-vec1) / (vec0 - vec1) / (vec4 - vec1));

                                   
                                    //{2,3,4,5}
                                    //determiant with basepoint 2 --> orientation "-"
                                    /*
                                    det(l1*(vec3 - vec0) + (vec0-vec2), l2*(vec4-vec1) + (vec1-vec2),
                                         l3*(vec5 - vec2))
                                    */

                                    Rational  lambda_1_lambda_2_lambda_3 =  det(vector2row(vec3-vec0) / (vec4 - vec1) / (vec5 - vec2));

                                    Rational  lambda_1_lambda_3 =  det(vector2row(vec3-vec0) / (vec1 - vec2) / (vec5 - vec2));
                                    Rational  lambda_2_lambda_3 =  det(vector2row(vec0-vec2) / (vec4 - vec1) / (vec5 - vec2));

                                    Rational  lambda_3 =  det(vector2row(vec0-vec2) / (vec1 - vec2) / (vec5 - vec2));

                                    if(det(vector2row(vec1-vec0) / (vec2-vec0) / (vec3-vec0)) > 0)
                                    {
                                        //1
                                        lobj_coeff[l1] +=  Rational(1,6)*lambda_1;


                                        //2
                                        qobj_coeff(l1,l2) += Rational(1,6)*lambda_1_lambda_2;
                                        lobj_coeff[l2] +=  Rational(1,6)*lambda_2;


                                        //3
                                        triobj_coeff[l1](l2, l3) +=  Rational(1,6)*lambda_1_lambda_2_lambda_3;
                                        qobj_coeff(l1,l3) += Rational(1,6)*lambda_1_lambda_3;
                                        qobj_coeff(l2,l3) += Rational(1,6)*lambda_2_lambda_3;
                                        lobj_coeff[l3] +=  Rational(1,6)*lambda_3;
      
                                                                                                 
                                    }
                                    else
                                    {
                                        //1
                                        lobj_coeff[l1] +=  Rational(1,6)*-lambda_1;


                                        //2
                                        qobj_coeff(l1,l2) += Rational(1,6)*-lambda_1_lambda_2;
                                        lobj_coeff[l2] +=  Rational(1,6)*-lambda_2;


                                        //3
                                        triobj_coeff[l1](l2, l3) +=  Rational(1,6)*-lambda_1_lambda_2_lambda_3;
                                        qobj_coeff(l1,l3) += Rational(1,6)*-lambda_1_lambda_3;
                                        qobj_coeff(l2,l3) += Rational(1,6)*-lambda_2_lambda_3;
                                        lobj_coeff[l3] +=  Rational(1,6)*-lambda_3;
                                                                                          
                                    }
                                    

                                    


                                }
                                        
                            }
                            else
                            {
                                Assert(constraints.size() == 0);

                                Vector<Rational> vec0 = all_vertices.row(vertices_so_far + triangle[0]);
                                Vector<Rational> vec1 = all_vertices.row(vertices_so_far + triangle[1]);
                                Vector<Rational> vec2 = all_vertices.row(vertices_so_far + triangle[2]);
                                Vector<Rational> vec3 = all_vertices.row(vertices_so_far + triangle[3]);
                                
                                

                                vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec3 = vec3.slice(range_from(ambient_dimension - slices_dimensions[j]));
                              

                                Rational Volume = cross_product3(vec1-vec0, vec2-vec0)*(vec3-vec0);

                                if(sub_points.contains(vertices_so_far + triangle[0]))
                                {
                                    FixVolume += Rational(1,6)*AbsoluteValue(Volume);
                                }
                                else
                                {
                                    ConterVolume += Rational(1,6)*AbsoluteValue(Volume);
                                }

                                

                                
                            }
                            

                        }

                        break;

                    case 4:
                        if(j == centerpoint_slice_index)
                        {
                            //TODO: 3d border with centerpoint
                            if(constraints.size() == 4)
                            {
                                //split vertices 2 by 2

                                /*

                                  1  ----   4
                                  | 2  --   | 5
                                  |/        |/
                                  0   ----  3

                                  NOTE: vec1 == vec4 and vec2 == vec5
                                */

                                int64 second_constraint_index = 0;

                                
                                for(int n = 1; n < constraints.size(); n++)
                                {
                                    if(constraints[0].start_vec_index != constraints[n].start_vec_index)
                                    {
                                        second_constraint_index = n;
                                    }

                                }


                                Vector<Rational> vec0 = all_vertices.row(constraints[0].start_vec_index);
                                Vector<Rational> vec1 = all_vertices.row(constraints[0].end_vec_index);
                               
                                int64 end_vec_index;
                                   
                                int64 add_edges = 0;
                               
                                for(int n = 1; n < constraints.size(); n++)
                                {
                                    if(constraints[0].start_vec_index == constraints[n].start_vec_index)
                                    {
                                        add_edges++;
                                        end_vec_index = n;
                                        break;
                                    }
                                }
                                Assert(add_edges == 1);

                                Vector<Rational> vec2 = all_vertices.row(end_vec_index);
                            
                                        
                                Vector<Rational> vec3 = all_vertices.row(constraints[second_constraint_index].start_vec_index);
                                Vector<Rational> vec4 = all_vertices.row(constraints[0].end_vec_index);
                                Vector<Rational> vec5 = all_vertices.row(end_vec_index);


                                //NOTE: THIS IS a slice form polymake!!
                                vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec3 = vec3.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec4 = vec4.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec5 = vec5.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                Vector<Rational> proj_centerpoint = centerpoint.slice(range_from(ambient_dimension - slices_dimensions[j]));

                                int64 l1 = edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index);
                                int64 l2 = edge_map(constraints[0].start_vec_index, end_vec_index);
                                int64 l4 = edge_map(constraints[second_constraint_index].start_vec_index, constraints[0].end_vec_index);

                                int64 l5 = edge_map(constraints[second_constraint_index].start_vec_index, end_vec_index);
                                        
                                        
                                //Triangulation
                                //{0,1,2,3},{1,2,3,4},{2,3,4,5}

                                //{0,1,2,3}
                                //determinant with base point 0 ---> orientation "-"
                                //det(l1*(vec1 - vec0), l2*(vec2-vec0), vec3 - vec0)
                                Rational lambda_1_lambda_2 =  det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec3 - vec0) / (proj_centerpoint - vec0));
                               

                                //{1,2,3,4}
                                //determinant with base point 3 ---> orientation "-"
                                //det(l1*(vec1 - vec0) - (vec3 - vec0), l2*(vec2-vec0)-
                                //    (vec3- vec0), l4(vec4 - vec 3))
                                        
                                       
                                Rational  lambda_1_lambda_2_lambda_4 =  det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec4 - vec3) / (proj_centerpoint - vec3));

                                //0: det(vec1 -vec 0, vec3-vec0, vec1-vec0) all in one plane!
                                Rational  lambda_1_lambda_4 =  det(vector2row(vec1-vec0) / -(vec3 - vec0) / (vec4 - vec3) / (proj_centerpoint - vec3));
                                        
                                Rational  lambda_2_lambda_4 =  det(vector2row(-(vec3-vec0)) / (vec2 - vec0) / (vec4 - vec3) / (proj_centerpoint - vec3));

                                //0:
                                Rational  lambda_4 =  det(vector2row(-(vec3-vec0)) / -(vec3 - vec0) / (vec3 - vec0) / (proj_centerpoint - vec3));

                                        


                                //{2,3,4,5}
                                //determiant with basepoint 3 --> orientation "+"
                                /*
                                  det(l2*(vec2 - vec0) - (vec3-vec0), l4*(vec4-vec3),
                                  l5*(vec5 - vec3))
                                */
                                Rational  lambda_2_lambda_4_lambda_5 =  det(vector2row(vec2-vec0) / (vec4 - vec3) / (vec5 - vec3) / (proj_centerpoint - vec3));


                                Rational  lambda_4_lambda_5 =  det(vector2row(-(vec3-vec0)) / (vec4 - vec3) / (vec5 - vec3) / (proj_centerpoint - vec3));

                                
                                //NOTE: Determin sign based on first triangle, this on is non singular since it arises from the triagle!
                                if(det(vector2row(vec1-vec0) / (vec2-vec0) / (vec3-vec0) / (proj_centerpoint - vec0)) > 0)  
                                {
                                    //1.
                                    qobj_coeff(l1,l2) += Rational(1,24)*lambda_1_lambda_2;
                                    

                                    //2.
                                    triobj_coeff[l1](l2, l4) +=  Rational(1,24)*lambda_1_lambda_2_lambda_4;
                                    qobj_coeff(l1,l4) += Rational(1,24)*lambda_1_lambda_4;

                                    qobj_coeff(l2,l4) += Rational(1,24)*lambda_2_lambda_4;

                                    lobj_coeff[l4] +=  Rational(1,24)*lambda_4;


                                    //3.
                                    triobj_coeff[l2](l4, l5) +=  Rational(1,24)*-lambda_2_lambda_4_lambda_5;
                                    qobj_coeff(l4,l5) += Rational(1,24)*-lambda_4_lambda_5;
                                                       
                                }
                                else
                                {
                                     //1.
                                    qobj_coeff(l1,l2) += Rational(1,24)*-lambda_1_lambda_2;
                                    

                                    //2.
                                    triobj_coeff[l1](l2, l4) +=  Rational(1,24)*-lambda_1_lambda_2_lambda_4;
                                    qobj_coeff(l1,l4) += Rational(1,24)*-lambda_1_lambda_4;

                                    qobj_coeff(l2,l4) += Rational(1,24)*-lambda_2_lambda_4;

                                    lobj_coeff[l4] +=  Rational(1,24)*-lambda_4;


                                    //3.
                                    triobj_coeff[l2](l4, l5) +=  Rational(1,24)*lambda_2_lambda_4_lambda_5;
                                    qobj_coeff(l4,l5) += Rational(1,24)*lambda_4_lambda_5;

                                }

                               

                                
                            }
                            else if(constraints.size() == 3)
                            {
                                if(constraints[0].start_vec_index == constraints[1].start_vec_index)
                                {
                                    /*

                                      1  
                                      | 2 
                                      |/        
                                      0--3

                                    */

                                    Vector<Rational> vec0 = all_vertices.row(constraints[0].start_vec_index);
                                    Vector<Rational> vec1 = all_vertices.row(constraints[0].end_vec_index);
                                    Vector<Rational> vec2 = all_vertices.row(constraints[1].end_vec_index);
                                    Vector<Rational> vec3 = all_vertices.row(constraints[2].end_vec_index);
                                    vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec3 = vec3.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    Vector<Rational> proj_centerpoint = centerpoint.slice(range_from(ambient_dimension - slices_dimensions[j]));

                                    int64 l1 = edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index);
                                    int64 l2 = edge_map(constraints[1].start_vec_index, constraints[1].end_vec_index);
                                    int64 l3 = edge_map(constraints[2].start_vec_index, constraints[2].end_vec_index);
                                    Rational lambda_1_lambda_2_lambda_3 =  det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec3 - vec0) / (proj_centerpoint - vec0));

                                    triobj_coeff[l1](l2, l3) +=  Rational(1,24)*AbsoluteValue(lambda_1_lambda_2_lambda_3);

                                }

                                else
                                {
                                     /*
                                       2 -
                                      | \  \
                                      |  |  \
                                      |  |   -
                                      |   1   \
                                      5   / \   -
                                      | 4    \   \
                                      |/       \  \
                                      x-----3-----0

                                    */

                                     //Triangulation
                                    //{0,1,2,3},{1,2,3,4},{2,3,4,5}
                                    //NOTE: WE have the same triangulation, but diffrent cutting edges

                                    Vector<Rational> vec0 = all_vertices.row(constraints[0].start_vec_index);
                                    Vector<Rational> vec1 = all_vertices.row(constraints[1].start_vec_index);
                                    Vector<Rational> vec2 = all_vertices.row(constraints[2].start_vec_index);
                                    Vector<Rational> vec3 = all_vertices.row(constraints[0].end_vec_index);
                                    Vector<Rational> vec4 = all_vertices.row(constraints[1].end_vec_index);
                                    Vector<Rational> vec5 = all_vertices.row(constraints[2].end_vec_index);

                                     //NOTE: vec3 == vec4 == vec5

                                    vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec3 = vec3.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec4 = vec4.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec5 = vec5.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    Vector<Rational> proj_centerpoint = centerpoint.slice(range_from(ambient_dimension - slices_dimensions[j]));

                                    int64 l1 = edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index);
                                    int64 l2 = edge_map(constraints[1].start_vec_index, constraints[1].end_vec_index);
                                    int64 l3 = edge_map(constraints[2].start_vec_index, constraints[2].end_vec_index);

                                    
                                    //{0,1,2,3}
                                    //determinant with base point 0 --> orientation "-"
                                    //det((vec1 - vec0), (vec2-vec0), l1*(vec3 - vec0))
                                    Rational lambda_1=  det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec3 - vec0) / (proj_centerpoint - vec0));

                                    //{1,2,3,4}
                                    //determinant with base point 1 --> orientation "-"
                                    //det((vec2 - vec1), l1*(vec3 - vec0) + (vec0 -vec1), l2*(vec4-vec1))
                                        
                                       
                                    Rational  lambda_1_lambda_2 =  det(vector2row(vec2-vec1) / (vec3 - vec0) / (vec4 - vec1) / (proj_centerpoint - vec1));

                                    Rational  lambda_2 =   det(vector2row(vec2-vec1) / (vec0 - vec1) / (vec4 - vec1) / (proj_centerpoint - vec1));

                                   
                                    //{2,3,4,5}
                                    //determiant with basepoint 2 --> orientation "-"
                                    /*
                                    det(l1*(vec3 - vec0) + (vec0-vec2), l2*(vec4-vec1) + (vec1-vec2),
                                         l3*(vec5 - vec2))
                                    */

                                    Rational  lambda_1_lambda_2_lambda_3 =  det(vector2row(vec3-vec0) / (vec4 - vec1) / (vec5 - vec2) / (proj_centerpoint - vec2));

                                    Rational  lambda_1_lambda_3 =  det(vector2row(vec3-vec0) / (vec1 - vec2) / (vec5 - vec2) / (proj_centerpoint - vec2));
                                    Rational  lambda_2_lambda_3 =  det(vector2row(vec0-vec2) / (vec4 - vec1) / (vec5 - vec2) / (proj_centerpoint - vec2));

                                    Rational  lambda_3 =  det(vector2row(vec0-vec2) / (vec1 - vec2) / (vec5 - vec2) / (proj_centerpoint - vec2));

                                    if(det(vector2row(vec1-vec0) / (vec2-vec0) / (vec3-vec0) / (proj_centerpoint - vec0)) > 0)
                                    {
                                        //1
                                        lobj_coeff[l1] +=  Rational(1,24)*lambda_1;


                                        //2
                                        qobj_coeff(l1,l2) += Rational(1,24)*lambda_1_lambda_2;
                                        lobj_coeff[l2] +=  Rational(1,24)*lambda_2;


                                        //3
                                        triobj_coeff[l1](l2, l3) +=  Rational(1,24)*lambda_1_lambda_2_lambda_3;
                                        qobj_coeff(l1,l3) += Rational(1,24)*lambda_1_lambda_3;
                                        qobj_coeff(l2,l3) += Rational(1,24)*lambda_2_lambda_3;
                                        lobj_coeff[l3] +=  Rational(1,24)*lambda_3;
      
                                                                                                 
                                    }
                                    else
                                    {
                                        //1
                                        lobj_coeff[l1] +=  Rational(1,24)*-lambda_1;


                                        //2
                                        qobj_coeff(l1,l2) += Rational(1,24)*-lambda_1_lambda_2;
                                        lobj_coeff[l2] +=  Rational(1,24)*-lambda_2;


                                        //3
                                        triobj_coeff[l1](l2, l3) +=  Rational(1,24)*-lambda_1_lambda_2_lambda_3;
                                        qobj_coeff(l1,l3) += Rational(1,24)*-lambda_1_lambda_3;
                                        qobj_coeff(l2,l3) += Rational(1,24)*-lambda_2_lambda_3;
                                        lobj_coeff[l3] +=  Rational(1,24)*-lambda_3;
                                                                                          
                                    }
                                    

                                    


                                }
                                        
                            }
                            else
                            {
                                Assert(constraints.size() == 0);

                                Vector<Rational> vec0 = all_vertices.row(vertices_so_far + triangle[0]);
                                Vector<Rational> vec1 = all_vertices.row(vertices_so_far + triangle[1]);
                                Vector<Rational> vec2 = all_vertices.row(vertices_so_far + triangle[2]);
                                Vector<Rational> vec3 = all_vertices.row(vertices_so_far + triangle[3]);
                                
                                

                                vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec3 = vec3.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                Vector<Rational> proj_centerpoint = centerpoint.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                 

                                Rational Volume = det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec3 - vec0) / (proj_centerpoint - vec0));;

                                if(sub_points.contains(vertices_so_far + triangle[0]))
                                {
                                    FixVolume += Rational(1,24)*AbsoluteValue(Volume);
                                }
                                else
                                {
                                    ConterVolume += Rational(1,24)*AbsoluteValue(Volume);
                                }

                                

                                
                            }

                        }
                        else
                        {

                            // We Always have 0, 4 or 6 contraits
                            if(constraints.size() == 6)
                            {
                                //2 / 3 or 3 / 2
                                int64 second_constraint_index = 0;
                                int64 third_constraint_index = 0;
                                

                                for(int n = 1; n < constraints.size(); n++)
                                {
                                    if((constraints[0].start_vec_index != constraints[n].start_vec_index)
                                       && (second_constraint_index == 0))
                                    {
                                        second_constraint_index = n;
                                    }

                                    else if(constraints[0].start_vec_index != constraints[n].start_vec_index)
                                    {
                                        third_constraint_index = n;
                                        
                                    }  

                                }

                                //Case 2 / 3
                                if(third_constraint_index == 0)
                                {
                                    Vector<Rational> vec0 = all_vertices.row(constraints[0].start_vec_index);
                                    Vector<Rational> vec1 = all_vertices.row(constraints[0].end_vec_index);

                                    int64 l1 = edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index);

                                    int64 end_vec_index[2];
                                   
                                    int64 add_edges = 0;

                                    for(int n = 1; n < constraints.size(); n++)
                                    {
                                        if(constraints[0].start_vec_index == constraints[n].start_vec_index)
                                        {
                                            add_edges++;
                                            end_vec_index[add_edges] = constraints[n].end_vec_index;
                                            
                                        }
                                        
                                    }
                                    Assert(add_edges == 2);

                                    Vector<Rational> vec2 = all_vertices.row(end_vec_index[0]);
                                    Vector<Rational> vec3 = all_vertices.row(end_vec_index[1]);

                                    int64 l2 = edge_map(constraints[0].start_vec_index, end_vec_index[0]);
                                    int64 l3 = edge_map(constraints[0].start_vec_index, end_vec_index[1]);
                                    

                                    Vector<Rational> vec4 = all_vertices.row(constraints[second_constraint_index].start_vec_index);

                                    //NOTE: vec5 == vec1 , vec6 == vec2 , vec7 == vec3
                                    Vector<Rational> vec5 = all_vertices.row(constraints[0].end_vec_index);

                                    int64 l5 = edge_map(constraints[second_constraint_index].start_vec_index, constraints[0].end_vec_index);

                                    Vector<Rational> vec6 = all_vertices.row(end_vec_index[0]);
                                    Vector<Rational> vec7 = all_vertices.row(end_vec_index[1]);

                                    int64 l6 = edge_map(constraints[second_constraint_index].start_vec_index, end_vec_index[0]);
                                    int64 l7 = edge_map(constraints[second_constraint_index].start_vec_index, end_vec_index[1]);



                                    //NOTE: THIS IS a slice form polymake!!
                                    vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec3 = vec3.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec4 = vec4.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec5 = vec5.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec6 = vec6.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec7 = vec7.slice(range_from(ambient_dimension - slices_dimensions[j]));

                                   
                                    //{{0 1 2 3 4} {1 2 3 4 5} {2 3 4 5 6} {3 4 5 6 7}}

                                    //{0,1,2,3,4}
                                    //determinant with base point 0 --> orientation "-"
                                    //det(l1*(vec1 - vec0), l2*(vec2-vec0), l3*(vec3 - vec0), (vec4 - vec0))
                                    Rational lambda_1_lambda_2_lambda_3 =  det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec3 - vec0) / (vec4 - vec0));

                                   
                                    //{1 2 3 4 5}
                                    //determinant with base point 4 --> orientation "-"
                                    //det(l1*(vec1 - vec0) + (vec0 - vec4) , l2*(vec2-vec0) + (vec0 - vec4),  l3*(vec3-vec0) + (vec0 - vec4), l5*(vec5 - vec4))
                                    Rational lambda_1_lambda_2_lambda_3_lambda_5 =  det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec3 - vec0) / (vec5 - vec4));
                                    
                                    Rational lambda_1_lambda_2_lambda_5 =  det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec0 - vec4) / (vec5 - vec4));
                                    Rational lambda_1_lambda_3_lambda_5 =  det(vector2row(vec1-vec0) / (vec0 - vec4) / (vec3 - vec0) / (vec5 - vec4));
                                    Rational lambda_2_lambda_3_lambda_5 =  det(vector2row(vec0-vec4) / (vec2 - vec0) / (vec3 - vec0) / (vec5 - vec4));

                                    //0:
                                    Rational lambda_1_lambda_5 =  det(vector2row(vec1-vec0) / (vec0 - vec4) / (vec0 - vec4) / (vec5 - vec4));
                                    Rational lambda_2_lambda_5 =  det(vector2row(vec0-vec4) / (vec2 - vec0) / (vec0 - vec4) / (vec5 - vec4));
                                    Rational lambda_3_lambda_5 =  det(vector2row(vec0-vec4) / (vec0 - vec4) / (vec3 - vec0) / (vec5 - vec4));

                                    //0:
                                    Rational lambda_5 =  det(vector2row(vec0-vec4) / (vec0 - vec4) / (vec0 - vec4) / (vec5 - vec4));
                                   

                                    // {2 3 4 5 6}
                                    //determinant with base point 4 --> orientation "-"
                                    //det(l2*(vec2 - vec0) + (vec0 - vec4) , l3*(vec3-vec0) + (vec0 - vec4),  l5*(vec5-vec4), l6*(vec6 - vec4))
                                    Rational lambda_2_lambda_3_lambda_5_lambda_6 =  det(vector2row(vec2-vec0) / (vec3 - vec0) / (vec5 - vec4) / (vec6 - vec4));
                                    
                                    Rational lambda_2_lambda_5_lambda_6 =  det(vector2row(vec2-vec0) / (vec0 - vec4) / (vec5 - vec4) / (vec6 - vec4));
                                    Rational lambda_3_lambda_5_lambda_6 =  det(vector2row(vec0-vec4) / (vec3 - vec0) / (vec5 - vec4) / (vec6 - vec4));
                                    
                                    Rational lambda_5_lambda_6 =  det(vector2row(vec0-vec4) / (vec0 - vec4) / (vec5 - vec4) / (vec6 - vec4));


                                    //{3 4 5 6 7}
                                    //determinant with base point 4 --> orientation "-"
                                    //det(l3*(vec3 - vec0) + (vec0 - vec4) , l5*(vec5-vec4), l6*(vec6-vec4), l7*(vec7 - vec4))
                                    Rational lambda_3_lambda_5_lambda_6_lambda_7 =  det(vector2row(vec3-vec0) / (vec5 - vec4) / (vec6 - vec4) / (vec7 - vec4));
                                    
                                    Rational lambda_5_lambda_6_lambda_7 =  det(vector2row(vec0-vec4) / (vec5 - vec4) / (vec6 - vec4) / (vec7 - vec4));
                                   
                                    if(det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec3 - vec0) / (vec4 - vec0)) > 0)
                                    {
                                        //1
                                        triobj_coeff[l1](l2, l3) +=  Rational(1,24)*lambda_1_lambda_2_lambda_3;

                                        //2
                                        quartobj_coeff[l1][l2](l3, l5) =  Rational(1,24)*lambda_1_lambda_2_lambda_3_lambda_5;
                                        triobj_coeff[l1](l2, l5) +=  Rational(1,24)*lambda_1_lambda_2_lambda_5;
                                        triobj_coeff[l1](l3, l5) +=  Rational(1,24)*lambda_1_lambda_3_lambda_5;
                                        triobj_coeff[l2](l3, l5) +=  Rational(1,24)*lambda_2_lambda_3_lambda_5;
                                        qobj_coeff(l1,l5) +=  Rational(1,24)*lambda_1_lambda_5;
                                        qobj_coeff(l2,l5) +=  Rational(1,24)*lambda_2_lambda_5;
                                        qobj_coeff(l3,l5) +=  Rational(1,24)*lambda_3_lambda_5;
                                        lobj_coeff[l5] +=  Rational(1,24)*lambda_5;

                                        //3
                                        quartobj_coeff[l2][l3](l5, l6) =  Rational(1,24)*lambda_2_lambda_3_lambda_5_lambda_6;
                                        triobj_coeff[l2](l5, l6) +=  Rational(1,24)*lambda_2_lambda_5_lambda_6;
                                        triobj_coeff[l3](l5, l6) +=  Rational(1,24)*lambda_3_lambda_5_lambda_6;
                                        qobj_coeff(l5,l6) +=  Rational(1,24)*lambda_5_lambda_6;
                                        
                                        //4
                                        quartobj_coeff[l3][l5](l6, l7) =  Rational(1,24)*lambda_3_lambda_5_lambda_6_lambda_7;
                                        triobj_coeff[l5](l6, l7) +=  Rational(1,24)*lambda_5_lambda_6_lambda_7;
                                        
                                       
                                    }
                                    else
                                    {
                                        //1
                                        triobj_coeff[l1](l2, l3) +=  Rational(1,24)*-lambda_1_lambda_2_lambda_3;

                                        //2
                                        quartobj_coeff[l1][l2](l3, l5) =  Rational(1,24)*-lambda_1_lambda_2_lambda_3_lambda_5;
                                        triobj_coeff[l1](l2, l5) +=  Rational(1,24)*-lambda_1_lambda_2_lambda_5;
                                        triobj_coeff[l1](l3, l5) +=  Rational(1,24)*-lambda_1_lambda_3_lambda_5;
                                        triobj_coeff[l2](l3, l5) +=  Rational(1,24)*-lambda_2_lambda_3_lambda_5;
                                        qobj_coeff(l1,l5) +=  Rational(1,24)*-lambda_1_lambda_5;
                                        qobj_coeff(l2,l5) +=  Rational(1,24)*-lambda_2_lambda_5;
                                        qobj_coeff(l3,l5) +=  Rational(1,24)*-lambda_3_lambda_5;
                                        lobj_coeff[l5] +=  Rational(1,24)*-lambda_5;

                                        //3
                                        quartobj_coeff[l2][l3](l5, l6) =  Rational(1,24)*-lambda_2_lambda_3_lambda_5_lambda_6;
                                        triobj_coeff[l2](l5, l6) +=  Rational(1,24)*-lambda_2_lambda_5_lambda_6;
                                        triobj_coeff[l3](l5, l6) +=  Rational(1,24)*-lambda_3_lambda_5_lambda_6;
                                        qobj_coeff(l5,l6) +=  Rational(1,24)*-lambda_5_lambda_6;
                                        
                                        //4
                                        quartobj_coeff[l3][l5](l6, l7) =  Rational(1,24)*-lambda_3_lambda_5_lambda_6_lambda_7;
                                        triobj_coeff[l5](l6, l7) +=  Rational(1,24)*-lambda_5_lambda_6_lambda_7;
                                        
                                    }

                                   

                                }
                                
                                //Case 3 / 2
                                
                                else
                                {
                                    Vector<Rational> vec0 = all_vertices.row(constraints[0].start_vec_index);
                                    Vector<Rational> vec1 = all_vertices.row(constraints[0].end_vec_index);

                                    int64 l1 = edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index);

                                    int64 end_vec_index;
                                   
                                    int64 add_edges = 0;

                                    for(int n = 1; n < constraints.size(); n++)
                                    {
                                        if(constraints[0].start_vec_index == constraints[n].start_vec_index)
                                        {
                                            add_edges++;
                                            end_vec_index = constraints[n].end_vec_index;
                                            
                                        }
                                        
                                    }
                                    Assert(add_edges == 1);

                                    Vector<Rational> vec2 = all_vertices.row(end_vec_index);

                                    int64 l2 = edge_map(constraints[0].start_vec_index, end_vec_index);

                                    


                                    Vector<Rational> vec3 = all_vertices.row(constraints[second_constraint_index].start_vec_index);
                                    Vector<Rational> vec4 = all_vertices.row(constraints[0].end_vec_index);

                                    int64 l4 = edge_map(constraints[second_constraint_index].start_vec_index, constraints[0].end_vec_index);


                                    Vector<Rational> vec5 = all_vertices.row(end_vec_index);

                                    int64 l5 = edge_map(constraints[second_constraint_index].start_vec_index, end_vec_index);

                                   

                                    Vector<Rational> vec6 = all_vertices.row(constraints[third_constraint_index].start_vec_index);
                                    Vector<Rational> vec7 = all_vertices.row(constraints[0].end_vec_index);

                                    int64 l7 = edge_map(constraints[third_constraint_index].start_vec_index, constraints[0].end_vec_index);


                                    Vector<Rational> vec8 = all_vertices.row(end_vec_index);

                                    int64 l8 = edge_map(constraints[third_constraint_index].start_vec_index, end_vec_index);

                                    
                                    // NOTE: vec4 == vec1 == vec7 , vec5 == vec2 == vec8 


                                    //NOTE: THIS IS a slice form polymake!!
                                    vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec3 = vec3.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec4 = vec4.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec5 = vec5.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec6 = vec6.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec7 = vec7.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec8 = vec8.slice(range_from(ambient_dimension - slices_dimensions[j]));

                                   
                                    //{{0 1 2 3 6} {1 2 3 4 6} {1 2 4 6 7} {2 3 4 5 6} {2 4 5 6 7} {2 5 6 7 8}}

                                    //{0,1,2,3,6}
                                    //determinant with base point 0 --> orientation "+"
                                    //det(l1*(vec1 - vec0), l2*(vec2-vec0), vec3 - vec0, vec6-vec0)

                                    Rational lambda_1_lambda_2 =  det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec3 - vec0) / (vec6 - vec0));
                                   


                                     //{1,2,3,4,6}
                                    //determinant with base point 3 --> orientation "+"
                                    //det(l1*(vec1 - vec0) + (vec0 - vec3), l2*(vec2-vec0) + (vec0-vec3), l4*(vec4 - vec3), vec6 - vec3)

                                    Rational lambda_1_lambda_2_lambda_4 =  det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec4 - vec3) / (vec6 - vec3));

                                    Rational lambda_1_lambda_4 =  det(vector2row(vec1-vec0) / (vec0 - vec3) / (vec4 - vec3) / (vec6 - vec3));

                                    Rational lambda_2_lambda_4 =  det(vector2row(vec0-vec3) / (vec2 - vec0) / (vec4 - vec3) / (vec6 - vec3));

                                    //0:
                                    Rational lambda_4 =  det(vector2row(vec0-vec3) / (vec0 - vec3) / (vec4 - vec3) / (vec6 - vec3));

 
                                    //{1,2,4,6,7}
                                    //determinant with base point 6 --> orientation "+"
                                    //det(l1*(vec1 - vec0) + (vec0 - vec6), l2*(vec2-vec0) + (vec0-vec6), l4*(vec4 - vec3) + (vec3 - vec6), l7*(vec7 - vec6))

                                    Rational lambda_1_lambda_2_lambda_4_lambda_7 =  det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec4 - vec3) / (vec7 - vec6));

                                    Rational lambda_1_lambda_2_lambda_7 =  det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec3 - vec6) / (vec7 - vec6));

                                    Rational lambda_1_lambda_4_lambda_7 =  det(vector2row(vec1-vec0) / (vec0 - vec6) / (vec4 - vec3) / (vec7 - vec6));

                                    Rational lambda_2_lambda_4_lambda_7 =  det(vector2row(vec0-vec6) / (vec2 - vec0) / (vec4 - vec3) / (vec7 - vec6));

                                    Rational lambda_1_lambda_7 =  det(vector2row(vec1-vec0) / (vec0 - vec6) / (vec3 - vec6) / (vec7 - vec6));

                                    Rational lambda_2_lambda_7 =  det(vector2row(vec0-vec6) / (vec2 - vec0) / (vec3 - vec6) / (vec7 - vec6));
                                    // 0:
                                    Rational lambda_4_lambda_7 =  det(vector2row(vec0-vec6) / (vec0 - vec6) / (vec4 - vec3) / (vec7 - vec6));

                                    //0:
                                    Rational lambda_7 =  det(vector2row(vec0-vec6) / (vec0 - vec6) / (vec3 - vec6) / (vec8 - vec6));


                                    // {2 3 4 5 6}
                                    //determinant with base point 3 --> orientation "-"
                                    //det(l2*(vec2 - vec0) + (vec0 - vec3), l4*(vec4-vec3), l5*(vec5 - vec3), (vec6 - vec3))

                                    Rational lambda_2_lambda_4_lambda_5 = det(vector2row(vec2-vec0) / (vec4 - vec3) / (vec5 - vec3) / (vec6 - vec3));
                                    //0:
                                    Rational lambda_4_lambda_5 =  det(vector2row(vec0-vec3) / (vec4 - vec3) / (vec5 - vec3) / (vec6 - vec3));

     
                                    //{2 4 5 6 7}
                                    //determinant with base point 6 --> orientation "-"
                                    //det(l2*(vec2 - vec0) + (vec0 - vec6), l4*(vec4 - vec3) + (vec3-vec6), l5*(vec5 - vec3) + (vec3 -vec6), l7*(vec7 - vec6))

                                    Rational lambda_2_lambda_4_lambda_5_lambda_7 = det(vector2row(vec2-vec0) / (vec4 - vec3) / (vec5 - vec3) / (vec7 - vec6));
                                    
                                    Rational two_lambda_2_lambda_4_lambda_7 =  det(vector2row(vec2-vec0) / (vec4 - vec3) / (vec3 - vec6) / (vec7 - vec6));

                                    Rational lambda_2_lambda_5_lambda_7 =  det(vector2row(vec2-vec0) / (vec3 - vec6) / (vec5 - vec3) / (vec7 - vec6));

                                    Rational lambda_4_lambda_5_lambda_7 =  det(vector2row(vec0-vec6) / (vec4 - vec3) / (vec5 - vec3) / (vec7 - vec6));

                                    Rational two_lambda_2_lambda_7 =  det(vector2row(vec2-vec0) / (vec3 - vec6) / (vec3 - vec6) / (vec7 - vec6));

                                    Rational two_lambda_4_lambda_7 =  det(vector2row(vec0-vec6) / (vec4 - vec6) / (vec3 - vec6) / (vec7 - vec6));

                                    Rational lambda_5_lambda_7 =  det(vector2row(vec0-vec6) / (vec3 - vec6) / (vec5 - vec3) / (vec7 - vec6));

                                    Rational two_lambda_7 =  det(vector2row(vec0-vec6) / (vec3 - vec6) / (vec3 - vec6) / (vec7 - vec6));

                                   
                                     //{2 5 6 7 8}
                                    //determinant with base point 6 --> orientation "+"
                                    //det(l2*(vec2 - vec0) + (vec0 - vec6), l5*(vec5 - vec3) + (vec3-vec6), l7*(vec7 - vec6), l8*(vec8 - vec6))

                                    Rational lambda_2_lambda_5_lambda_7_lambda_8 = det(vector2row(vec2-vec0) / (vec5 - vec3) / (vec7 - vec6) / (vec8 - vec6));
                                    
                                    Rational lambda_2_lambda_7_lambda_8 =  det(vector2row(vec2-vec0) / (vec3 - vec6) / (vec7 - vec6) / (vec8 - vec6));

                                    Rational lambda_5_lambda_7_lambda_8 =  det(vector2row(vec0-vec6) / (vec5 - vec3) / (vec7 - vec3) / (vec8 - vec6));

                                    Rational lambda_7_lambda_8 =  det(vector2row(vec0-vec6) / (vec3 - vec6) / (vec7 - vec6) / (vec8 - vec6));


                                   
                                    if(det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec3 - vec0) / (vec6 - vec0)) > 0)
                                    {
                                        //1
                                         qobj_coeff(l1,l2) +=  Rational(1,24)*lambda_1_lambda_2;
                                        //2
                                        triobj_coeff[l1](l2, l4) +=  Rational(1,24)*lambda_1_lambda_2_lambda_4;
                                        qobj_coeff(l1,l4) +=  Rational(1,24)*lambda_1_lambda_4;
                                        qobj_coeff(l2,l4) +=  Rational(1,24)*lambda_2_lambda_4;
                                        lobj_coeff[l4] +=  Rational(1,24)*lambda_4;

                                        //3
                                        quartobj_coeff[l1][l2](l4, l7) =  Rational(1,24)*lambda_1_lambda_2_lambda_4_lambda_7;
                                        triobj_coeff[l1](l2, l7) +=  Rational(1,24)*lambda_1_lambda_2_lambda_7;
                                        triobj_coeff[l1](l4, l7) +=  Rational(1,24)*lambda_1_lambda_4_lambda_7;
                                        triobj_coeff[l2](l4, l7) +=  Rational(1,24)*lambda_2_lambda_4_lambda_7;
                                        qobj_coeff(l1,l7) +=  Rational(1,24)*lambda_1_lambda_7;
                                        qobj_coeff(l2,l7) +=  Rational(1,24)*lambda_2_lambda_7;
                                        qobj_coeff(l2,l7) +=  Rational(1,24)*lambda_4_lambda_7;
                                        lobj_coeff[l7] +=  Rational(1,24)*lambda_7;

                                        //4
                                        triobj_coeff[l2](l4, l5) +=  Rational(1,24)*-lambda_2_lambda_4_lambda_5;
                                        qobj_coeff(l4,l5) +=  Rational(1,24)*-lambda_4_lambda_5;

                                        //5

                                        quartobj_coeff[l2][l4](l5, l7) =  Rational(1,24)*-lambda_2_lambda_4_lambda_5_lambda_7;
                                        triobj_coeff[l2](l4, l7) +=  Rational(1,24)*-two_lambda_2_lambda_4_lambda_7;
                                        triobj_coeff[l2](l5, l7) +=  Rational(1,24)*-lambda_2_lambda_4_lambda_7;
                                        triobj_coeff[l4](l5, l7) +=  Rational(1,24)*-lambda_2_lambda_4_lambda_7;
                                        qobj_coeff(l2,l7) +=  Rational(1,24)*-two_lambda_2_lambda_7;
                                        qobj_coeff(l4,l7) +=  Rational(1,24)*-two_lambda_4_lambda_7;
                                        qobj_coeff(l5,l7) +=  Rational(1,24)*-lambda_5_lambda_7;
                                        lobj_coeff[l7] +=  Rational(1,24)*-two_lambda_7;
                                        
                                        //6
                                        quartobj_coeff[l2][l5](l7, l8) =  Rational(1,24)*lambda_2_lambda_5_lambda_7_lambda_8;
                                        triobj_coeff[l2](l7, l8) +=  Rational(1,24)*lambda_2_lambda_7_lambda_8;
                                        triobj_coeff[l5](l7, l8) +=  Rational(1,24)*lambda_5_lambda_7_lambda_8;
                                        qobj_coeff(l7,l8) +=  Rational(1,24)*lambda_7_lambda_8;
                                           
                                    }

                                    else
                                    {                                       
                                        //1
                                         qobj_coeff(l1,l2) +=  Rational(1,24)*-lambda_1_lambda_2;
                                        //2
                                        triobj_coeff[l1](l2, l4) +=  Rational(1,24)*-lambda_1_lambda_2_lambda_4;
                                        qobj_coeff(l1,l4) +=  Rational(1,24)*-lambda_1_lambda_4;
                                        qobj_coeff(l2,l4) +=  Rational(1,24)*-lambda_2_lambda_4;
                                        lobj_coeff[l4] +=  Rational(1,24)*-lambda_4;

                                        //3
                                        quartobj_coeff[l1][l2](l4, l7) =  Rational(1,24)*-lambda_1_lambda_2_lambda_4_lambda_7;
                                        triobj_coeff[l1](l2, l7) +=  Rational(1,24)*-lambda_1_lambda_2_lambda_7;
                                        triobj_coeff[l1](l4, l7) +=  Rational(1,24)*-lambda_1_lambda_4_lambda_7;
                                        triobj_coeff[l2](l4, l7) +=  Rational(1,24)*-lambda_2_lambda_4_lambda_7;
                                        qobj_coeff(l1,l7) +=  Rational(1,24)*-lambda_1_lambda_7;
                                        qobj_coeff(l2,l7) +=  Rational(1,24)*-lambda_2_lambda_7;
                                        qobj_coeff(l2,l7) +=  Rational(1,24)*-lambda_4_lambda_7;
                                        lobj_coeff[l7] +=  Rational(1,24)*-lambda_7;

                                        //4
                                        triobj_coeff[l2](l4, l5) +=  Rational(1,24)*lambda_2_lambda_4_lambda_5;
                                        qobj_coeff(l4,l5) +=  Rational(1,24)*lambda_4_lambda_5;

                                        //5

                                        quartobj_coeff[l2][l4](l5, l7) =  Rational(1,24)*lambda_2_lambda_4_lambda_5_lambda_7;
                                        triobj_coeff[l2](l4, l7) +=  Rational(1,24)*two_lambda_2_lambda_4_lambda_7;
                                        triobj_coeff[l2](l5, l7) +=  Rational(1,24)*lambda_2_lambda_4_lambda_7;
                                        triobj_coeff[l4](l5, l7) +=  Rational(1,24)*lambda_2_lambda_4_lambda_7;
                                        qobj_coeff(l2,l7) +=  Rational(1,24)*two_lambda_2_lambda_7;
                                        qobj_coeff(l4,l7) +=  Rational(1,24)*two_lambda_4_lambda_7;
                                        qobj_coeff(l5,l7) +=  Rational(1,24)*lambda_5_lambda_7;
                                        lobj_coeff[l7] +=  Rational(1,24)*two_lambda_7;
                                        
                                        //6
                                        quartobj_coeff[l2][l5](l7, l8) =  Rational(1,24)*-lambda_2_lambda_5_lambda_7_lambda_8;
                                        triobj_coeff[l2](l7, l8) +=  Rational(1,24)*-lambda_2_lambda_7_lambda_8;
                                        triobj_coeff[l5](l7, l8) +=  Rational(1,24)*-lambda_5_lambda_7_lambda_8;
                                        qobj_coeff(l7,l8) +=  Rational(1,24)*-lambda_7_lambda_8;

                                    }
      
                                }

                            }

                            else if(constraints.size() == 4)
                            {
                                //1 / 4 or 4 / 1

                                //Case 1 / 4
                                if(constraints[0].start_vec_index == constraints[1].start_vec_index)
                                {
                                    Vector<Rational> vec0 = all_vertices.row(constraints[0].start_vec_index);
                                    Vector<Rational> vec1 = all_vertices.row(constraints[0].end_vec_index);

                                    int64 l1 = edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index);

                                    Vector<Rational> vec2 = all_vertices.row(constraints[1].end_vec_index);

                                    int64 l2 = edge_map(constraints[0].start_vec_index, constraints[1].end_vec_index);

                                    Vector<Rational> vec3 = all_vertices.row(constraints[2].end_vec_index);

                                    int64 l3 = edge_map(constraints[0].start_vec_index, constraints[2].end_vec_index);

                                    Vector<Rational> vec4 = all_vertices.row(constraints[3].end_vec_index);

                                    int64 l4 = edge_map(constraints[0].start_vec_index, constraints[3].end_vec_index);

                                    vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec3 = vec3.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec4 = vec4.slice(range_from(ambient_dimension - slices_dimensions[j]));

                                    //{0,1,2,3,4}
                                    //determinant with base point 0
                                    //det(l1*(vec1 - vec0), l2*(vec2-vec0), l3*(vec3 - vec0), l4*(vec4-vec0))

                                    Rational lambda_1_lambda_2_lambda_3_lamdba_4 =  det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec3 - vec0) / (vec4 - vec0));

                                    quartobj_coeff[l1][l2](l3, l4) +=  Rational(1,24)*AbsoluteValue(lambda_1_lambda_2_lambda_3_lamdba_4);


                                }

                                //Case 4 / 1
                                else
                                {
                                    Vector<Rational> vec0 = all_vertices.row(constraints[0].start_vec_index);
                                    Vector<Rational> vec1 = all_vertices.row(constraints[0].end_vec_index);

                                    int64 l1 = edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index);

                                    Vector<Rational> vec2 = all_vertices.row(constraints[1].start_vec_index);
                                    Vector<Rational> vec3 = all_vertices.row(constraints[1].end_vec_index);

                                    int64 l2 = edge_map(constraints[1].start_vec_index, constraints[1].end_vec_index);

                                    Vector<Rational> vec4 = all_vertices.row(constraints[2].start_vec_index);
                                    Vector<Rational> vec5 = all_vertices.row(constraints[2].end_vec_index);

                                    int64 l3 = edge_map(constraints[2].start_vec_index, constraints[2].end_vec_index);

                                    Vector<Rational> vec6 = all_vertices.row(constraints[3].start_vec_index);
                                    Vector<Rational> vec7 = all_vertices.row(constraints[3].end_vec_index);

                                    int64 l4 = edge_map(constraints[3].start_vec_index, constraints[3].end_vec_index);

                                    //NOTE: vec1 == vec3 == vec5 == vec7
                                    
                                    vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec3 = vec3.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec4 = vec4.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec5 = vec5.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec6 = vec6.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec7 = vec7.slice(range_from(ambient_dimension - slices_dimensions[j]));

                                    //{{0 1 2 4 6} {1 2 3 4 6} {1 3 4 5 6} {1 3 5 6 7}}


                                    //{0 1 2 4 6}
                                    //determinant with base point 0 --> orientation "-"
                                    //det(l1*(vec1 - vec0), (vec2-vec0), (vec4 - vec0), (vec6-vec0))

                                    Rational lambda_1 =  det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec4 - vec0) / (vec6 - vec0));


                                    //{1 2 3 4 6}
                                    //determinant with base point 2 --> orientation "-"
                                    //det(l1*(vec1 - vec0) + (vec0 - vec2), l2*(vec3-vec2), (vec4 - vec2), (vec6-vec2))

                                    Rational lambda_1_lambda_2 =  det(vector2row(vec1-vec0) / (vec3 - vec2) / (vec4 - vec2) / (vec6 - vec2));
                                    Rational lambda_2 =  det(vector2row(vec0-vec2) / (vec3 - vec2) / (vec4 - vec2) / (vec6 - vec2));
    
                                   
                                    //{1 3 4 5 6}
                                    //determinant with base point --> orientation "-"
                                    //det(l1*(vec1 - vec0) + (vec0 - vec4), l2*(vec3-vec2) + (vec2 - vec4), l3*(vec5 - vec4), l4*(vec7-vec6) + (vec6 - vec4))

                                    Rational lambda_1_lambda_2_lambda_3_lambda_4 =  det(vector2row(vec1-vec0) / (vec3 - vec2) / (vec5 - vec4) / (vec7 - vec6));
                                    Rational lambda_1_lambda_2_lambda_3 =  det(vector2row(vec1-vec0) / (vec3 - vec2) / (vec5 - vec4) / (vec6 - vec4));
                                    Rational lambda_1_lambda_3_lambda_4 =  det(vector2row(vec1-vec0) / (vec2 - vec4) / (vec5 - vec4) / (vec7 - vec6));
                                    Rational lambda_2_lambda_3_lambda_4 =  det(vector2row(vec0-vec4) / (vec3 - vec2) / (vec5 - vec4) / (vec7 - vec6));

                                    Rational lambda_1_lambda_3=  det(vector2row(vec1-vec0) / (vec2 - vec4) / (vec5 - vec4) / (vec6 - vec4));
                                    Rational lambda_2_lambda_3 =  det(vector2row(vec0-vec4) / (vec3 - vec2) / (vec5 - vec4) / (vec6 - vec4));
                                    Rational lambda_3_lambda_4 =  det(vector2row(vec0-vec4) / (vec2 - vec4) / (vec5 - vec4) / (vec7 - vec6));

                                    Rational lambda_3 =  det(vector2row(vec0-vec4) / (vec2 - vec4) / (vec5 - vec4) / (vec6 - vec4));
                                   


                                    //{1 3 5 6 7}
                                    //determinant with base point 6 --> orientation "-"
                                    //det(l1*(vec1 - vec0) + (vec0 - vec6), l2*(vec3-vec2) + (vec2 - vec6), l3*(vec5 - vec4) + (vec4 - vec6), l4*(vec7-vec6))

                                    Rational two_lambda_1_lambda_2_lambda_3_lambda_4 =  det(vector2row(vec1-vec0) / (vec3 - vec2) / (vec5 - vec4) / (vec7 - vec6));
                                    Rational lambda_1_lambda_2_lambda_4 =  det(vector2row(vec1-vec0) / (vec3 - vec2) / (vec4 - vec6) / (vec6 - vec4));
                                    Rational two_lambda_1_lambda_3_lambda_4 =  det(vector2row(vec1-vec0) / (vec3 - vec6) / (vec5 - vec4) / (vec7 - vec6));
                                    Rational two_lambda_2_lambda_3_lambda_4 =  det(vector2row(vec0-vec6) / (vec3 - vec2) / (vec5 - vec4) / (vec7 - vec6));

                                    Rational lambda_1_lambda_4=  det(vector2row(vec1-vec0) / (vec2 - vec6) / (vec4 - vec6) / (vec7 - vec6));
                                    Rational lambda_2_lambda_4 =  det(vector2row(vec0-vec6) / (vec3 - vec2) / (vec4 - vec6) / (vec7 - vec6));
                                    Rational two_lambda_3_lambda_4 =  det(vector2row(vec0-vec6) / (vec2 - vec6) / (vec5 - vec4) / (vec7 - vec6));

                                    Rational lambda_4 =  det(vector2row(vec0-vec4) / (vec2 - vec4) / (vec5 - vec4) / (vec6 - vec4));
                                     
                                    if(det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec4 - vec0) / (vec6 - vec0)) > 0)
                                    {
                                        //1
                                        lobj_coeff[l1] +=  Rational(1,24)*lambda_1;

                                        //2
                                        qobj_coeff(l1,l2) +=  Rational(1,24)*lambda_1_lambda_2;
                                        lobj_coeff[l2] +=  Rational(1,24)*lambda_2;

                                        //3
                                        quartobj_coeff[l1][l2](l3, l4) =  Rational(1,24)*lambda_1_lambda_2_lambda_3_lambda_4;
                                        triobj_coeff[l1](l2, l3) +=  Rational(1,24)*lambda_1_lambda_2_lambda_3;
                                        triobj_coeff[l1](l3, l4) +=  Rational(1,24)*lambda_1_lambda_3_lambda_4;
                                        triobj_coeff[l2](l3, l4) +=  Rational(1,24)*lambda_2_lambda_3_lambda_4;
                                        qobj_coeff(l1,l3) +=  Rational(1,24)*lambda_1_lambda_3;
                                        qobj_coeff(l2,l3) +=  Rational(1,24)*lambda_2_lambda_3;
                                        qobj_coeff(l3,l4) +=  Rational(1,24)*lambda_3_lambda_4;
                                        lobj_coeff[l3] +=  Rational(1,24)*lambda_3;
                                        
                                        //4
                                        quartobj_coeff[l1][l2](l3, l4) =  Rational(1,24)*two_lambda_1_lambda_2_lambda_3_lambda_4;
                                        triobj_coeff[l1](l2, l4) +=  Rational(1,24)*lambda_1_lambda_2_lambda_4;
                                        triobj_coeff[l1](l3, l4) +=  Rational(1,24)*two_lambda_1_lambda_3_lambda_4;
                                        triobj_coeff[l2](l3, l4) +=  Rational(1,24)*two_lambda_2_lambda_3_lambda_4;
                                        qobj_coeff(l1,l4) +=  Rational(1,24)*lambda_1_lambda_4;
                                        qobj_coeff(l2,l4) +=  Rational(1,24)*lambda_2_lambda_4;
                                        qobj_coeff(l3,l4) +=  Rational(1,24)*two_lambda_3_lambda_4;
                                        lobj_coeff[l4] +=  Rational(1,24)*lambda_4;
                                    }
                                    else
                                    {
                                         //1
                                        lobj_coeff[l1] +=  Rational(1,24)*-lambda_1;

                                        //2
                                        qobj_coeff(l1,l2) +=  Rational(1,24)*-lambda_1_lambda_2;
                                        lobj_coeff[l2] +=  Rational(1,24)*-lambda_2;

                                        //3
                                        quartobj_coeff[l1][l2](l3, l4) =  Rational(1,24)*-lambda_1_lambda_2_lambda_3_lambda_4;
                                        triobj_coeff[l1](l2, l3) +=  Rational(1,24)*-lambda_1_lambda_2_lambda_3;
                                        triobj_coeff[l1](l3, l4) +=  Rational(1,24)*-lambda_1_lambda_3_lambda_4;
                                        triobj_coeff[l2](l3, l4) +=  Rational(1,24)*-lambda_2_lambda_3_lambda_4;
                                        qobj_coeff(l1,l3) +=  Rational(1,24)*-lambda_1_lambda_3;
                                        qobj_coeff(l2,l3) +=  Rational(1,24)*-lambda_2_lambda_3;
                                        qobj_coeff(l3,l4) +=  Rational(1,24)*-lambda_3_lambda_4;
                                        lobj_coeff[l3] +=  Rational(1,24)*-lambda_3;
                                        
                                        //4
                                        quartobj_coeff[l1][l2](l3, l4) =  Rational(1,24)*-two_lambda_1_lambda_2_lambda_3_lambda_4;
                                        triobj_coeff[l1](l2, l4) +=  Rational(1,24)*-lambda_1_lambda_2_lambda_4;
                                        triobj_coeff[l1](l3, l4) +=  Rational(1,24)*-two_lambda_1_lambda_3_lambda_4;
                                        triobj_coeff[l2](l3, l4) +=  Rational(1,24)*-two_lambda_2_lambda_3_lambda_4;
                                        qobj_coeff(l1,l4) +=  Rational(1,24)*-lambda_1_lambda_4;
                                        qobj_coeff(l2,l4) +=  Rational(1,24)*-lambda_2_lambda_4;
                                        qobj_coeff(l3,l4) +=  Rational(1,24)*-two_lambda_3_lambda_4;
                                        lobj_coeff[l4] +=  Rational(1,24)*-lambda_4;
                                    }

                                }

                            }

                            else
                            {
                                Assert(constraints.size() == 0);

                                Vector<Rational> vec0 = all_vertices.row(vertices_so_far + triangle[0]);
                                Vector<Rational> vec1 = all_vertices.row(vertices_so_far + triangle[1]);
                                Vector<Rational> vec2 = all_vertices.row(vertices_so_far + triangle[2]);
                                Vector<Rational> vec3 = all_vertices.row(vertices_so_far + triangle[3]);
                                Vector<Rational> vec4 = all_vertices.row(vertices_so_far + triangle[4]);
                                
                                

                                vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec3 = vec3.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                vec4 = vec4.slice(range_from(ambient_dimension - slices_dimensions[j]));
                              

                                Rational Volume = det(vector2row(vec1-vec0) / (vec2 - vec0) / (vec3 - vec0) / (vec4 - vec0));

                                if(sub_points.contains(vertices_so_far + triangle[0]))
                                {
                                    FixVolume += Rational(1,24)*AbsoluteValue(Volume);
                                }
                                else
                                {
                                    ConterVolume += Rational(1,24)*AbsoluteValue(Volume);
                                }

                            }

                        }
                        
                        break;
                }
                

            }

            vertices_so_far += slices_vertices[j].rows();
        }


        

        //Build optimization problem



        
        bool32 clamp_precision = false;
        double    sol[number_lambda + ambient_dimension];
        double    objval;
      
          

#if OPT_DEBUG            
        cout << "The model objective is: " << endl;
        cout << convert_to<double>(qobj_coeff) << endl;
        cout << convert_to<double>(lobj_coeff) << endl;
        cout << convert_to<double>(cobj_coeff) << endl;
#endif


REPEAT:       //gurobi test:

        GRBenv   *env   = NULL;
        GRBmodel *model = NULL;
        int       error = 0;
        double    obj_lin[number_lambda];
        int       qrow[ambient_dimension];
        int       qcol[ambient_dimension];
        double    qval[ambient_dimension];
        int       lind[ambient_dimension];
        double    lval[ambient_dimension];
        int       optimstatus;
       
        double    lambda_lb[number_lambda];
        double    lambda_ub[number_lambda];
        double    u_lb[ambient_dimension];
        double    u_ub[ambient_dimension];
        int       u2_qrow[ambient_dimension];
        int       u2_qcol[ambient_dimension];
        double    u2_qval[ambient_dimension];

        //NOTE: The bounds on the lambda are very important otherwise the computation takes forever... 
              for(int32 w = 0; w < number_lambda; w++)
          {
              lambda_lb[w] = 0;
              lambda_ub[w] = 1;
          }

           
        for(int32 w = 0; w < ambient_dimension; w++)
        {
            u_lb[w] = -1+1e-3;
            u_ub[w] = 1-1e-3;
        }
           
            
         

        for(int32 w = 0; w < lobj_coeff.size(); w++)
        {
            obj_lin[w] =  convert_to<double>(lobj_coeff[w]);
        }

            
            
#if SOL_GUROBI                                      

        /* Create environment */

        error = GRBloadenv(&env, NULL);
        if (error) goto QUIT1;

#if OPT_NO_OUTPUT            
        // Stop gurobi from printing to the console
        error = GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0);
        if (error) goto QUIT1;
#endif


        error = GRBsetdblparam(env, GRB_DBL_PAR_TIMELIMIT, 10.0);
        if (error) goto QUIT1;

        /* Create an empty model */

        error = GRBnewmodel(env, &model, "qcp", 0, NULL, NULL, NULL, NULL, NULL);
        if (error) goto QUIT1;


        /* Add variables */
            
        // Lambda + linear objective function term
        // 
            error = GRBaddvars(model, number_lambda, 0, NULL, NULL, NULL, obj_lin, lambda_lb, lambda_ub, NULL , NULL);
        
            if (error) goto QUIT1;

            // u
            //
            error = GRBaddvars(model, ambient_dimension, 0, NULL, NULL, NULL, NULL, u_lb, u_ub, NULL, NULL);
            if (error) goto QUIT1;


            /* Add constraints */

            // Hyperplane constraint for each lambda
            //
            for(int32 k = 0; k < edge_map.rows(); k++)
            {
                for(int32 l = 0; l < edge_map.cols(); l++)
                {
                   
                    
                    if(edge_map(k,l) != -1)
                    {
                        Vector<Rational> vec0 = all_vertices.row(k);
                        Vector<Rational> vec1 = all_vertices.row(l);
                        Vector<Rational> num_vec = centerpoint - vec0;
                        Vector<Rational> denum_vec = vec1 - vec0;

                        
                        for(int32 v = 0; v < ambient_dimension; v++)
                        {
                            //Note number_lambda + v is the u indicex
                            qrow[v] = edge_map(k,l);
                            qcol[v] = number_lambda + v;
                            qval[v] = convert_to<double>(denum_vec[v]);

                            lind[v] = number_lambda+v;
                            lval[v] = - convert_to<double>(num_vec[v]);

                            
                           
                        }

#if OPT_DEBUG
                
                        cout << "lambd for " << k << " " << l <<endl;
                        cout << "num =" <<  convert_to<double>(num_vec) << endl;
                        cout << "denum =" <<  convert_to<double>(denum_vec) << endl;
#endif
                
                        error = GRBaddqconstr(model, ambient_dimension, lind, lval, ambient_dimension, qrow, qcol, qval, GRB_EQUAL, 0, NULL);
                        if (error) goto QUIT1;
                    }
                }

            }
            
            
               

            // Cone constraint for u
            //
            for(int32 w = 0; w < cone_h_representation.rows(); w++)
            {

                for(int32 v = 0; v < ambient_dimension; v++)
                {
                    //Note number_lambda + v is the u index
                    if(clamp_precision == true)
                    {
                        lind[v] = number_lambda+v;
                        lval[v] = decimal_truncate(convert_to<double>(cone_h_representation(w,v)),3);
                        
                    }
                    else
                    {
                        lind[v] = number_lambda+v; lval[v] = convert_to<double>(cone_h_representation(w,v));
                    }
                
                }
                

                error = GRBaddconstr(model, ambient_dimension, lind, lval, GRB_GREATER_EQUAL, 0.0, NULL);
                if (error) goto QUIT1;

            }



            
            //Norm constraint for u, or just u > 0
            //
            for(int32 w = 0; w < ambient_dimension; w++)
            {
                //Note number_lambda + w is the u index
                u2_qrow[w] = number_lambda+w; u2_qcol[w] = number_lambda+w; u2_qval[w]= 1.0;
            }
            
            error = GRBaddqconstr(model, 0, NULL, NULL, ambient_dimension, u2_qrow, u2_qcol, u2_qval, GRB_EQUAL, 1.0, "norm");
            if (error) goto QUIT1;
            



            
            /* mimization */

            error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE);
            if (error) goto QUIT1;

            
            //Add quadratic objective terms

            
            for (int32 w = 0; w < qobj_coeff.rows(); w++)
            {
                //We need for number_lambda = 4:
                //int obj_qrow[number_lambda] = {w,w,w,w};
                //int obj_qcol[number_lambda] = {0,1,2,3};
                int obj_qrow[number_lambda];
                int obj_qcol[number_lambda];
                double obj_qval[number_lambda];

                for(int32 v = 0; v < qobj_coeff.cols(); v++)
                {
                    obj_qrow[v] = w;
                    obj_qcol[v] = v;
                    obj_qval[v] = convert_to<double>(qobj_coeff(w,v));
                }
       

                error = GRBaddqpterms(model, number_lambda, obj_qrow, obj_qcol, obj_qval);
                if (error) goto QUIT1;
            }
            
            
            //THIS IS DEBUG TEST:
            if((iteration == 4) && (i == 2))
            {
                cout << "HERE We want to STOP!" << endl;
                error = GRBwrite(model, "debug_qcp_21_06.lp");
                
                if (error)
                {
                    cout << "Problem when writing error model" << endl;
                };

                //Assert(1 == 0);
            }

            

            /* Optimize model */

            error = GRBoptimize(model);
            if (error) goto QUIT1;

            /* Capture solution information */

            error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
            if (error) goto QUIT1;


            

            printf("\nOptimization complete\n");
            if (optimstatus == GRB_OPTIMAL)
            {
                error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
                if (error) goto QUIT1;

                error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, number_lambda+ambient_dimension, sol);
                if (error) goto QUIT1;
                
                printf("Optimal objective: %.4e\n", objval);

                for(int32 w = 0; w < number_lambda; w++)
                {
                    cout << "l" << w << "=" << sol[w]  << " ";
                   
                }

                //NOTE: WE cast the direction back to Rational; Do I have to worry about precision??
                
                for(int32 w = 0; w < ambient_dimension; w++)
                {                   
                    minimal_directions(i,w) = convert_to<Rational>(sol[number_lambda + w]);
                    cout << "u" << w << "=" << sol[number_lambda + w]  << " ";
                }

                cout << std::defaultfloat << std::setprecision(6);
                cout << endl;

                /* Free model */

                GRBfreemodel(model);

                /* Free environment */

                GRBfreeenv(env);

                
            }
            else if (optimstatus == GRB_INF_OR_UNBD)
            {
                gurobi_infeasibilities++;
                //Assert(1 == 0);

                objval = convert_to<double>(TotalVolume);
                for(int w = 0; w < number_lambda+ambient_dimension; w++)
                {
                    sol[w] = 0.0;
                }
                printf("Model is infeasible or unbounded\n");
                
            }
            else
            {
                gurobi_infeasibilities++;

                //Assert(1 == 0); 

                objval = convert_to<double>(TotalVolume);
                for(int w = 0; w < number_lambda+ambient_dimension; w++)
                {
                    sol[w] = 0.0;
                    }
                printf("Optimization was stopped early\n");
            }

            
        
QUIT1:

            /* Error reporting */

            if (error)
            {
                printf("ERROR: %s\n", GRBgeterrormsg(env));
                
                error = GRBwrite(model, "error_qcp.lp");
                if (error)
                {
                    cout << "Problem when writing error model" << endl;
                };

                exit(1);
               
            }

            /* Free model */

            GRBfreemodel(model);

            /* Free environment */

            GRBfreeenv(env);

           

            
#endif
            
#if SOL_BARON
            //BARON FILE:


            //NOTE: This creates or overwrites the file
            std::ofstream file("opt_model.bar"); 

            if (!file)
            {
                std::cerr << "Error: could not open file for writing.\n";
                Assert(1 == 0);
            }
            else
            {
                file << "POSITIVE_VARIABLES  ";
                for (int w = 0; w < number_lambda; w++)
                {
                    file << "x" << w;
                    if (w < number_lambda - 1)
                    {
                        file << ",";
                    }
                }

                    

                file << ";\n";

                file << "\n";

                file << "VARIABLES  ";
                for (int w = number_lambda; w < number_lambda + ambient_dimension; w++)
                {
                    file << "x" << w;
                    if (w < number_lambda + ambient_dimension - 1)
                    {
                        file << ",";
                    }
                }

                file << ";\n";

                file << "\n"; 

                // Write LOWER_BOUNDS 
                       file << "LOWER_BOUNDS{\n";
                for (int w = 0; w < number_lambda; w++)
                {
                    file << "x" << w << ": " << std::setprecision(16) << lambda_lb[w] << ";\n";
                        
                
                }
                for (int w = 0; w < ambient_dimension; w++)
                {
                    /*
                      double x = u_lb[w];

                      std::stringstream ss;
                      //ss << std::setprecision( std::numeric_limits<double>::digits10+2);
                      ss << std::setprecision( std::numeric_limits<int>::max() );
                      ss << x;

                      std::string str = ss.str();
                    */
                    file << "x" << (number_lambda + w) << ": " << std::setprecision(16) << u_lb[w] << ";\n";
                
                        
                }
                file << "}\n";

                file << "\n"; 

                // Write UPPER_BOUNDS block
                file << "UPPER_BOUNDS{\n";
                for (int w = 0; w < number_lambda; w++)
                {
                    file << "x" << w << ": " << std::setprecision(16) << lambda_ub[w] << ";\n";
                }
                for (int w = 0; w < ambient_dimension; w++)
                {
                    file << "x" <<  (number_lambda + w)  << ": " << std::setprecision(16) << u_ub[w] << ";\n";
                }
                file << "}\n";

                    
                file << "\n";

                file << "EQUATIONS  ";
                for (int w = 0; w < number_lambda + cone_h_representation.rows() + 1; w++)
                {
                    file << "e" << w;
                    if (w < number_lambda + cone_h_representation.rows() + 1 - 1)
                    {
                        file << ",";
                    }
                }
                file << ";\n";

                file << "\n";
                file << "\n";
                    

                // Hyperplane constraint for each lambda
                //
                int num_lambda_equations = 0;
                for(int32 k = 0; k < edge_map.rows(); k++)
                {
                    for(int32 l = 0; l < edge_map.cols(); l++)
                    {
                   
                    
                        if(edge_map(k,l) != -1)
                        {
                            file << "e" << num_lambda_equations << ": ";
                            num_lambda_equations++;
                            Vector<Rational> vec0 = all_vertices.row(k);
                            Vector<Rational> vec1 = all_vertices.row(l);
                            Vector<Rational> num_vec = centerpoint - vec0;
                            Vector<Rational> denum_vec = vec1 - vec0;


                            //quadratic
                            for(int32 v = 0; v < ambient_dimension; v++)
                            {
                                //Note number_lambda + v is the u indicex

                                double coeff = convert_to<double>(denum_vec[v]);
                                file << std::showpos;
                                 
                                file << std::setprecision(16) << coeff;
 
                                file << std::noshowpos;
                                file << "*x" <<  edge_map(k,l) << "*x" << (number_lambda + v) << " ";
                                  
                            
                           
                            }

                            for(int32 v = 0; v < ambient_dimension; v++)
                            {
                                //Note number_lambda + v is the u indicex

                                double coeff = - convert_to<double>(num_vec[v]);
                                file  << std::showpos;
                                 
                                file  << std::setprecision(16) << coeff;
 
                                file  << std::noshowpos;
                                file <<  "*x" << (number_lambda + v) << " ";
                           
                            }

                            file << "== 0.0;" << endl;
                            file << "\n" << endl;

                        }
                    }

                }

                Assert(number_lambda == num_lambda_equations);

                // Cone constraint for u
                //
                int num_cone_equations = num_lambda_equations;
                for(int32 w = 0; w < cone_h_representation.rows(); w++)
                {

                    file << "e" << num_cone_equations << ": ";
                    num_cone_equations++;
                    for(int32 v = 0; v < ambient_dimension; v++)
                    {
                        double coeff = convert_to<double>(cone_h_representation(w,v));
                        file << std::showpos;
                                 
                        file << std::setprecision(16) << coeff;
 
                        file << std::noshowpos;
                        file << "*x" << (number_lambda + v) << " ";
                            
                
                    }

                    file << ">= 0.0;" << endl;
                    file << "\n" << endl;

                }
                    

                //Norm constraint for u, or just u > 0
                //
                file << "e" << number_lambda +  cone_h_representation.rows() << ": ";
                for(int32 w = 0; w < ambient_dimension; w++)
                {

                    double coeff = 1.0;
                    file << std::showpos;
                                 
                    file << std::setprecision(16) << coeff;
 
                    file << std::noshowpos;
                    file << "*x" << (number_lambda + w) << "*x" << (number_lambda + w) << " ";
                }
            
                    
                file << "== 1.0;" << endl;
                file << "\n" << endl;

                    
                file << "OBJ: minimize    ";

                    
                //quadratic objective:
                for (int32 w = 0; w < qobj_coeff.rows(); w++)
                {
                    for(int32 v = 0; v < qobj_coeff.cols(); v++)
                    {
                            
                        double coeff = convert_to<double>(qobj_coeff(w,v));
                        if(coeff != 0)
                        {
                            file << std::showpos;
                                 
                            file << std::setprecision(16) << coeff;
 
                            file << std::noshowpos;
                            file << "*x" <<  w << "*x" << v << " ";
                                
                        }
                            
                    }
       
                }

                //linear objective:
                for (int32 w = 0; w < lobj_coeff.size(); w++)
                {
                    double coeff = convert_to<double>(lobj_coeff[w]);
                    if(coeff != 0)
                    {
                        file << std::showpos;
                                 
                        file << std::setprecision(16) << coeff;
 
                        file << std::noshowpos;
                        file << "*x" <<  w;
                        if(w != (lobj_coeff.size()-1))
                        {
                            file <<  " ";
                        }
                                
                    }
 
                }
                    
                file << ";" << endl;
                file.close();


            }

            system("../baron-lin64/baron opt_model.bar"); 

            read_baron_results(&objval, sol, number_lambda + ambient_dimension);

            //NOTE: WE cast the direction back to Rational; Do I have to worry about precision??
                
            for(int32 w = 0; w < ambient_dimension; w++)
            {                   
                minimal_directions(i,w) = convert_to<Rational>(sol[number_lambda + w]);
                cout << "u" << w << "=" << sol[number_lambda + w]  << " ";
            }

#endif
            
#if 0
            
            vertices_so_far = 0;
           

            Matrix<double> real_vertices = convert_to<double>(all_vertices);

            for(int j = 0; j < num_slices; j++)
            {
                Array<Set<Int>> slice_triangulation = slices_triangulations[j];

            


                for(int k = 0; k < slice_triangulation.size(); k++)
                {

                    std::vector<lambda_u_constraint> constraints;
                    Vector<Int> triangle = Vector<Int>(slice_triangulation[k]);
                    //compute the total volume of the triangle:


                    for(int l = 0; l < triangle.size(); l++)
                    {
                        if(sub_boundary_points.contains(vertices_so_far + triangle[l]))
                        {
                            for(int m = 0; m < triangle.size(); m++)
                            {
                                if(conter_boundary_points.contains(vertices_so_far + triangle[m]))
                                {
                                    lambda_u_constraint cutting_edge = {};
                                    cutting_edge.start_vec_index = vertices_so_far + triangle[l];
                                    cutting_edge.end_vec_index = vertices_so_far + triangle[m];
                                    constraints.push_back(cutting_edge);
                           
                                }

                            }
                       
                        }
                
                    }

                
                    //NOTE: HERE THE TRIANGULATION IS A BOUNDARY TRIANGULATION

                    switch (slices_dimensions[j])
                    {
                        case 1:
                            // We Always have 0 or 1 contraits
                            if(j == centerpoint_slice_index)
                            {
                                if(constraints.size() > 0)
                                {
                                    //NOTHING TO DO HERE
                                }
                                else
                                {
                                    //This case should never happen!! -- The triangulation should leave a
                                    // line as it is
                                    //NOTE: MAYBE if the centerpoint is the vertex??
                                    Assert(1 == 0);
                                }
                            }
                            else
                            {
                                if(constraints.size() > 0)
                                {
                                    Assert(constraints.size() == 1);
                                    Vector<double> vec0 = real_vertices.row(constraints[0].start_vec_index);
                                    Vector<double> vec1 = real_vertices.row(constraints[0].end_vec_index);

                                    

                                    double lambda1 = sol[edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index)];

                                    Vector<double> intersection_vec_1 = real_vertices.row(constraints[0].start_vec_index) + lambda1*(real_vertices.row(constraints[0].end_vec_index) - real_vertices.row(constraints[0].start_vec_index));
                                    //NOTE: THIS IS a slice form polymake!!
                                    vec0 = vec0.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    intersection_vec_1 = intersection_vec_1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    Assert(vec0.dim() == 1);

                                    ConterVolume += AbsoluteValue(vec1[0] - intersection_vec_1[0]);
                                }
                                else
                                {
                                    //NOTHING TO DO HERE
                                }
                            
                            }
                            break;
                            
                        case 2:
                            // We Always have 0 or 1 contraits

                            if(j == centerpoint_slice_index)
                            {
                                if(constraints.size() > 0)
                                {
                                    Assert(constraints.size() == 1);

                                    Vector<double> vec1 = real_vertices.row(constraints[0].end_vec_index);

                                    
                                    double lambda1 = sol[edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index)];
                                  
                                    Vector<double> intersection_vec_1 = real_vertices.row(constraints[0].start_vec_index) + lambda1*(real_vertices.row(constraints[0].end_vec_index) - real_vertices.row(constraints[0].start_vec_index));
                                    
                               
                                   

                                    //NOTE: THIS IS a slice form polymake!!
                                    vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    intersection_vec_1 = intersection_vec_1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    Vector<Rational> proj_centerpoint = centerpoint.slice(range_from(ambient_dimension - slices_dimensions[j]));

                                    Matrix<double> real_conter_vertices = vector2row(vec1);
                                    if(lambda1 != 1)
                                    {
                                        real_conter_vertices /= intersection_vec_1;
                                    }
                                    real_conter_vertices /= convert_to<double>(proj_centerpoint);
                                          
                                    
                                    Matrix<Rational> rational_conter_vertices = convert_to<Rational>(real_conter_vertices);

                                    BigObject conter_triangle("Polytope<Rational>");

                                    conter_triangle.take("VERTICES") << (ones_vector<Rational>() | rational_conter_vertices);
             
                                    Rational ConterVolume_triangle = conter_triangle.give("VOLUME");

                                    ConterVolume += ConterVolume_triangle;

                                }
                                else
                                {

                                    //NOTHING TO DO
                                
                                }
                            }
                            else
                            {
                                if(constraints.size() > 0)
                                {
                                    Assert(constraints.size() == 2);
                                    
                                    Vector<double> vec1 = real_vertices.row(constraints[0].end_vec_index);
                                    Vector<double> vec2 = real_vertices.row(constraints[1].end_vec_index);
                                    

                                    double lambda1 = sol[edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index)];
                                    double lambda2 = sol[edge_map(constraints[1].start_vec_index, constraints[1].end_vec_index)];

                                    Vector<double> intersection_vec_1 = real_vertices.row(constraints[0].start_vec_index) + lambda1*(real_vertices.row(constraints[0].end_vec_index) - real_vertices.row(constraints[0].start_vec_index));
                                    Vector<double> intersection_vec_2 = real_vertices.row(constraints[1].start_vec_index) + lambda2*(real_vertices.row(constraints[1].end_vec_index) - real_vertices.row(constraints[1].start_vec_index));

                                    vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    intersection_vec_1 = intersection_vec_1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    intersection_vec_2 = intersection_vec_2.slice(range_from(ambient_dimension - slices_dimensions[j]));

                                    //NOTE: MAYBE it is better to avoid doing this case destinction and
                                    //just do "POINTS" instead of "VERTICES"
                                    if(constraints[0].start_vec_index == constraints[1].start_vec_index)
                                    {
                                        /* case:
                                           current |----------| current_conter
                                           | \        |
                                           |   \      |
                                           |      \   |
                                           |        \ |
                                           next |----------| next_conter

                                        */


                                        Matrix<double> real_conter_vertices = vector2row(intersection_vec_1);
                                        if((lambda1 != 0) || (lambda2 != 0))
                                        {
                                            real_conter_vertices /= intersection_vec_2;
                                        }
                                        if(lambda1 != 1)
                                        {
                                            real_conter_vertices /= vec1;
                                        }
                                        if(lambda2 != 1)
                                        {
                                            real_conter_vertices /= vec2;
                                        }

                                        Matrix<Rational> rational_conter_vertices = convert_to<Rational>(real_conter_vertices);

                                        BigObject conter_triangle("Polytope<Rational>");

                                        conter_triangle.take("VERTICES") << (ones_vector<Rational>() | rational_conter_vertices);
             
                                        Rational ConterVolume_triangle = conter_triangle.give("VOLUME");

                                        ConterVolume += ConterVolume_triangle;

                                        
                               
                                

                                    }

                                    else
                                    {

                                        /* case:
                                           current |----------| current_conter
                                           |        / |
                                           |      /   |
                                           |   /      |
                                           | /        |
                                           next |----------| next_conter

                                        */

                                        Matrix<double> real_conter_vertices = vector2row(intersection_vec_1);
                                        if((lambda1 != 1) || (lambda2 != 1))
                                        {
                                            real_conter_vertices /= intersection_vec_2;
                                            real_conter_vertices /= vec1;
                                        }
                                        
                                    
                                        Matrix<Rational> rational_conter_vertices = convert_to<Rational>(real_conter_vertices);

                                        BigObject conter_triangle("Polytope<Rational>");

                                        conter_triangle.take("VERTICES") << (ones_vector<Rational>() | rational_conter_vertices);
             
                                        Rational ConterVolume_triangle = conter_triangle.give("VOLUME");

                                        ConterVolume += ConterVolume_triangle;

                                    }
                                }
                                else
                                {

                                    //NOTHING TO DO
                                
                                }

                            }
                            break;
                        case 3:

                            if(j == centerpoint_slice_index)
                            {
                                // We Always have 0 or 2 contraits
                                if(constraints.size() > 0)
                                {
                                    Assert(constraints.size() == 2);
                                    
                                    Vector<double> vec1 = real_vertices.row(constraints[0].end_vec_index);
                                    Vector<double> vec2 = real_vertices.row(constraints[1].end_vec_index);
                                    

                                    double lambda1 = sol[edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index)];
                                    double lambda2 = sol[edge_map(constraints[1].start_vec_index, constraints[1].end_vec_index)];

                                    Vector<double> intersection_vec_1 = real_vertices.row(constraints[0].start_vec_index) + lambda1*(real_vertices.row(constraints[0].end_vec_index) - real_vertices.row(constraints[0].start_vec_index));
                                    Vector<double> intersection_vec_2 = real_vertices.row(constraints[1].start_vec_index) + lambda2*(real_vertices.row(constraints[1].end_vec_index) - real_vertices.row(constraints[1].start_vec_index));

                                    vec1 = vec1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    vec2 = vec2.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    intersection_vec_1 = intersection_vec_1.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    intersection_vec_2 = intersection_vec_2.slice(range_from(ambient_dimension - slices_dimensions[j]));
                                    Vector<Rational> proj_centerpoint = centerpoint.slice(range_from(ambient_dimension - slices_dimensions[j]));

                                    //NOTE: MAYBE it is better to avoid doing this case destinction and
                                    //just do "POINTS" instead of "VERTICES"
                                    if(constraints[0].start_vec_index == constraints[1].start_vec_index)
                                    {
                                        /* case:
                                           current |----------| current_conter
                                           | \        |
                                           |   \      |
                                           |      \   |
                                           |        \ |
                                           next |----------| next_conter

                                        */


                                        Matrix<double> real_conter_vertices = vector2row(intersection_vec_1);
                                        if((lambda1 != 0) || (lambda2 != 0))
                                        {
                                            real_conter_vertices /= intersection_vec_2;
                                        }
                                        if(lambda1 != 1)
                                        {
                                            real_conter_vertices /= vec1;
                                        }
                                        if(lambda2 != 1)
                                        {
                                            real_conter_vertices /= vec2;
                                        }

                                        real_conter_vertices /= convert_to<double>(proj_centerpoint);

                                        Matrix<Rational> rational_conter_vertices = convert_to<Rational>(real_conter_vertices);

                                        BigObject conter_triangle("Polytope<Rational>");

                                        conter_triangle.take("VERTICES") << (ones_vector<Rational>() | rational_conter_vertices);
             
                                        Rational ConterVolume_triangle = conter_triangle.give("VOLUME");

                                        ConterVolume += ConterVolume_triangle;

                                        
                               
                                

                                    }

                                    else
                                    {

                                        /* case:
                                           current |----------| current_conter
                                           |        / |
                                           |      /   |
                                           |   /      |
                                           | /        |
                                           next |----------| next_conter

                                        */

                                        Matrix<double> real_conter_vertices = vector2row(intersection_vec_1);
                                        if((lambda1 != 1) || (lambda2 != 1))
                                        {
                                            real_conter_vertices /= intersection_vec_2;
                                            real_conter_vertices /= vec1;
                                        }

                                        real_conter_vertices /= convert_to<double>(proj_centerpoint);
                                    
                                        Matrix<Rational> rational_conter_vertices = convert_to<Rational>(real_conter_vertices);

                                        BigObject conter_triangle("Polytope<Rational>");

                                        conter_triangle.take("VERTICES") << (ones_vector<Rational>() | rational_conter_vertices);
             
                                        Rational ConterVolume_triangle = conter_triangle.give("VOLUME");

                                        ConterVolume += ConterVolume_triangle;
                                    }
                    
                                }

                                else
                                {
                                    //NOTHING TO DO

                                }
                            }
                            else
                            {
                                // We Always have 0, 3 or 4 contraits
                                if(constraints.size() == 4)
                                {
                                    //TODO
                                }
                            

                            }

                            break;
                            
                            

                            
                    }

                }
                vertices_so_far += slices_vertices[j].rows();
            }

            
#endif
           


            
            


           

            
            real64 sub_volume = objval + convert_to<double>(FixVolume);

            cout << "number_lambda: " << number_lambda << endl;
            cout << "ConterVolume: " << convert_to<double>(ConterVolume) << endl;
           
            minimal_volume_signature_indices.push_back(i);
            minimal_volumes.push_back(sub_volume);
           
            
            real64 total_volume = convert_to<double>(ConterVolume) + objval + convert_to<double>(FixVolume);
            
            cout << "Volume: " << convert_to<double>(ConterVolume) << "(conter) + " << objval << "(objval) + " << convert_to<double>(FixVolume) << "(fix-sub) = " << total_volume << endl;

            cout << "Total Volume " << convert_to<double>(TotalVolume) << endl;


            /*
            if((total_volume > 167) || (total_volume < 166))
            {
                cout << "THE COMPUTED VOLUMES ARE NOT CORRECT !!!!";
                Assert(1 == 0);
                   
            }
            */
        
       

        

       
       




/*
        
        BigObject PointConfig("PointConfiguration<Rational>");

        Matrix<Rational> points = sub_boundary_vertices / conter_boundary_vertices / centerpoint;

        PointConfig.take("POINTS") << (ones_vector<Rational>() | points);

        //topaz::GeometricSimplicialComplex<Rational> triangulation = PointConfig.give("TRIANGULATION");

        Array<Set<Set<Int>>> triangulations= call_function("topcom_all_triangulations", PointConfig);


        cout << triangulations[0] << endl;

        
*/





        

        /* 
        if((total_volume > 1001) || (total_volume < 999))
        {
            cout << "THE COMPUTED VOLUMES ARE NOT CORRECT !!!!";
            Assert(1 == 0);
                   
        }
        */
         
           
    } 

    


    real64  minimal_value = minimal_volumes[0];
    int64 minimal_value_index = 0;
    for(int i = 0; i < minimal_volume_signature_indices.size(); i++)
    {

        Set<Int> sig_points;
        
        for(int j = 0; j < signatures.cols(); j++)
        {
            
            if(signatures(minimal_volume_signature_indices[i],j))
            {
                sig_points.collect(j);   
            }
          
        }
        
        cout << sig_points << "  Volume:" << minimal_volumes[i] <<  endl;

        if(minimal_volumes[i] < minimal_value)
        {
            minimal_value = minimal_volumes[i];
            minimal_value_index = i;
        }
    }

    Set<Int> minimal_sub_points;
        
    for(int j = 0; j < signatures.cols(); j++)
    {
            
        if(signatures(minimal_volume_signature_indices[minimal_value_index],j))
        {
            minimal_sub_points.collect(j);   
        }
          
    }

    cout << endl;
    
    cout << "MINIMAL VOLUME " << minimal_value << " attaind at " << minimal_sub_points << endl;

    halfspace Halfspace = {};
    Halfspace.direction = minimal_directions.row(minimal_value_index);
    Halfspace.value = minimal_volumes[minimal_value_index];
    Halfspace.infeasibilities = gurobi_infeasibilities;

    cout << "Minimal Direction: " << convert_to<double>(Halfspace.direction) << endl;

    return Halfspace;
}







log_point sub_grad_descend_slices(Array<Matrix<Rational>> slices_vertices, Vector<Rational> point,  int slice_index, memory_arena *Arena)
{

    int centerpoint_slice_index = slice_index;

    int64 ambient_dimension = slices_vertices[0].cols();


    int num_slices = slices_vertices.size();

    Matrix<Rational> all_vertices;
    Array<Array<Set<Int>>> slices_triangulations(num_slices);
    Array<Array<Set<Int>>> slices_edges(num_slices);

    Rational TotalVolume = 0;
    Array<Rational> slices_volumes(num_slices);
    Array<Int> slices_dimensions(num_slices);

     
    Matrix<Rational> centerpoint_slice_Inequalities;
    int64 centerpoint_slice_dimension;

    for(int i = 0; i < num_slices; i++)
    {


        all_vertices /= slices_vertices[i];
        BigObject slice("Polytope<Rational>");
        slice.take("VERTICES") << (ones_vector<Rational>() | slices_vertices[i]);

        Matrix<Rational> Inequalities = slice.give("FACETS");

        int64 dimension = slice.call_method("DIM");
        slices_dimensions[i] = dimension;


        
        //Build the inequality system (NOTE: The first vector corresponds to the right-hand-side)
        

        if(ambient_dimension - dimension == 1)
        {
            Vector<Rational> affine_eq_x_1 = zero_vector<Rational>(ambient_dimension + 1);
            Vector<Rational> affine_eq_x_2 = zero_vector<Rational>(ambient_dimension  + 1);
            affine_eq_x_1[1] = -1;
            affine_eq_x_2[1] = 1;
            affine_eq_x_1[0] = slices_vertices[i](0,0);
            affine_eq_x_2[0] = -slices_vertices[i](0,0);
            Inequalities /= affine_eq_x_1;
            Inequalities /= affine_eq_x_2;

            
            /*
              Vector<Rational> affine_eq_1 = {0,-1, 0, 0};
              Vector<Rational> affine_eq_2 = {0, 1, 0, 0};
              affine_eq_1[0] = slices_vertices[i](0,0);
              affine_eq_2[0] = -slices_vertices[i](0,0);
              Inequalities /= affine_eq_1;
              Inequalities /= affine_eq_2;
            */
        }
        else if(ambient_dimension - dimension == 2)
        {
            Vector<Rational> affine_eq_x_1 = zero_vector<Rational>(ambient_dimension + 1);
            Vector<Rational> affine_eq_x_2 = zero_vector<Rational>(ambient_dimension + 1);
            affine_eq_x_1[1] = -1;
            affine_eq_x_2[1] = 1;
            affine_eq_x_1[0] = slices_vertices[i](0,0);
            affine_eq_x_2[0] = -slices_vertices[i](0,0);
            Inequalities /= affine_eq_x_1;
            Inequalities /= affine_eq_x_2;

            Vector<Rational> affine_eq_y_1 = zero_vector<Rational>(ambient_dimension + 1);
            Vector<Rational> affine_eq_y_2 = zero_vector<Rational>(ambient_dimension + 1);
            affine_eq_y_1[2] = -1;
            affine_eq_y_2[2] = 1;
            affine_eq_y_1[0] = slices_vertices[i](0,1);
            affine_eq_y_2[0] = -slices_vertices[i](0,1);
            Inequalities /= affine_eq_y_1;
            Inequalities /= affine_eq_y_2;
        }

        else if(ambient_dimension - dimension == 3)
        {
            Vector<Rational> affine_eq_x_1 = zero_vector<Rational>(ambient_dimension + 1);
            Vector<Rational> affine_eq_x_2 = zero_vector<Rational>(ambient_dimension + 1);
            affine_eq_x_1[1] = -1;
            affine_eq_x_2[1] = 1;
            affine_eq_x_1[0] = slices_vertices[i](0,0);
            affine_eq_x_2[0] = -slices_vertices[i](0,0);
            Inequalities /= affine_eq_x_1;
            Inequalities /= affine_eq_x_2;

            Vector<Rational> affine_eq_y_1 = zero_vector<Rational>(ambient_dimension + 1);
            Vector<Rational> affine_eq_y_2 = zero_vector<Rational>(ambient_dimension + 1);
            affine_eq_y_1[2] = -1;
            affine_eq_y_2[2] = 1;
            affine_eq_y_1[0] = slices_vertices[i](0,1);
            affine_eq_y_2[0] = -slices_vertices[i](0,1);
            Inequalities /= affine_eq_y_1;
            Inequalities /= affine_eq_y_2;

            Vector<Rational> affine_eq_z_1 = zero_vector<Rational>(ambient_dimension + 1);
            Vector<Rational> affine_eq_z_2 = zero_vector<Rational>(ambient_dimension + 1);
            affine_eq_z_1[3] = -1;
            affine_eq_z_2[3] = 1;
            affine_eq_z_1[0] = slices_vertices[i](0,2);
            affine_eq_z_2[0] = -slices_vertices[i](0,2);
            Inequalities /= affine_eq_z_1;
            Inequalities /= affine_eq_z_2;

        }
        else if(ambient_dimension - dimension > 3)
        {
            cout << "The case that n(integer) > 3 currently isn't handled" << endl;
            Assert(false);

            //Debugger doesn't like that for some reason
            //Matrix<Rational> slice_affine_h_rep = slice.give("AFFINE_HULL");
        
            //Inequalities /= slice_affine_h_rep;
        }
   

        cout << Inequalities << endl;

        if(i == centerpoint_slice_index)
        {
            centerpoint_slice_Inequalities = Inequalities;
            centerpoint_slice_dimension = dimension;
        
            if(!is_feasible(&Inequalities, point))
            {
                cout << "Given point " << point << " not in slice " << centerpoint_slice_index << endl;
                Assert(1 == 0);
            }
        }

        if((slices_vertices[i].cols() - dimension > 0) && (dimension > 0))
        {
            //NOTE: This methode is save to call (i. e. we don't get any irrational data)
            // Since the slices are axis aligned!
            Map<Rational,Rational> rel_volume = slice.give("RELATIVE_VOLUME");
            Assert(rel_volume.size() == 1);
            slices_volumes[i] = rel_volume[1];
            TotalVolume += rel_volume[1];    
        }
        else
        {
            Rational Volume = slice.give("VOLUME");
            slices_volumes[i] = Volume;
            TotalVolume += Volume;
        }

        
        if(is_feasible(&Inequalities, point))
        {
            cout << "FEASABLE for Slice: " << i << endl;
            //cout << "affine_h_rep" << slice_affine_h_rep << endl;
            centerpoint_slice_index = i;
            BigObject simplicalcomplex;

            slice.give("TRIANGULATION") >> simplicalcomplex;

            BigObject simplicalcomplex_boundary;

            simplicalcomplex.give("BOUNDARY") >> simplicalcomplex_boundary;

            Array<Set<Int>> boundary_triangulation = simplicalcomplex_boundary.give("FACETS");
            slices_triangulations[i] = boundary_triangulation;

            BigObject triangulation_graph("Graph<Undirected>");
            triangulation_graph = simplicalcomplex_boundary.give("GRAPH");
       
            Array<Set<Int>> triangulation_edges = triangulation_graph.call_method("EDGES");

            slices_edges[i] = triangulation_edges;
        }
        else
        {
            BigObject simplicalcomplex;

            slice.give("TRIANGULATION") >> simplicalcomplex;

            Array<Set<Int>> triangulation = simplicalcomplex.give("FACETS");
            slices_triangulations[i] = triangulation;

            BigObject triangulation_graph("Graph<Undirected>");
            triangulation_graph = simplicalcomplex.give("GRAPH");
       
            Array<Set<Int>> triangulation_edges = triangulation_graph.call_method("EDGES");

            slices_edges[i] = triangulation_edges;

            
        }
    }


     // End of preprocessing data!..

    
    int64 Num_Iteration = 30;
    Vector<Rational> current_point = point;
    halfspace start_halfspace = minimum_halfspace_slices(slices_vertices, all_vertices, slices_triangulations, slices_edges, TotalVolume, slices_volumes, slices_dimensions, ambient_dimension, centerpoint_slice_index, current_point, 0);
    Vector<Rational> current_direction = - start_halfspace.direction;
    
    //Rescale the direction in the current slice to be unit norm

    Rational scale_denum = 1;

    for(int i = 0; i < ambient_dimension - centerpoint_slice_dimension; i++)
    {
        scale_denum -= current_direction[i]*current_direction[i];
    }

    Assert(scale_denum >= 0);

    if(scale_denum > 0)
    {
        for(int j = 0; j < current_direction.size(); j++)
        {
            //Maybe its better to leave the direction in double?
            if(j >= ambient_dimension - centerpoint_slice_dimension)
            {
                double coordinate = convert_to<double>(current_direction[j]);
                double new_coordinate = coordinate / std::sqrt(convert_to<double>(scale_denum));
                current_direction[j] = convert_to<Rational>(new_coordinate);
            }
            else
            {
                current_direction[j] = 0;
            }
            
                
        }
    }

    
    int64 unfeasable_count = 0;
    
   

    log_point *first_point = PushStruct(Arena, log_point);
    first_point->point.X = convert_to<real64>(current_point[0]);
    first_point->point.Y = convert_to<real64>(current_point[1]);
    first_point->point.Z = convert_to<real64>(current_point[2]);
    first_point->direction.X = convert_to<real64>(current_direction[0]);
    first_point->direction.Y = convert_to<real64>(current_direction[1]);
    first_point->direction.Z = convert_to<real64>(current_direction[2]);
    first_point->value = start_halfspace.value;
    first_point->total_infeasibilities = start_halfspace.infeasibilities;

    
    for (int i = 1; i < Num_Iteration + 1; i++)
    {
        Vector<Rational> test_point = current_point + Rational(1,i)*current_direction;

        cout << convert_to<double>(test_point) << endl;

        //Check for feasible
        //NOTE: The direction vector is noramlized!
        int64 outer_loop_iteration = 0;

        if(!is_feasible(&centerpoint_slice_Inequalities, test_point))
        {
            unfeasable_count++;

            int64 accuracy = 1;
            Rational lambda =  Rational(1, accuracy);
      
            while(!is_feasible(&centerpoint_slice_Inequalities, test_point))
            {
                cout << convert_to<double>(test_point) << endl;
                accuracy++;
                lambda = Rational(1, accuracy);
                test_point = current_point + lambda*Rational(1,i)*current_direction;
            }


            Assert(is_feasible(&centerpoint_slice_Inequalities, test_point));

        }

        /*
        //NOTE: This procedure was implemented before the scaling of the cone-h-representation
        // The new code is in the beginning of the minimum_halfspace_slices function
        // It helped in some cases... but has no logical reason to work
          Vector<double> d_test_point = convert_to<double>(test_point);
                
          Vector<Rational> rounded_test_point = zero_vector<Rational>(ambient_dimension);

        
          rounded_test_point[0] = convert_to<Rational>(round_to_digits(d_test_point[0],9));
          rounded_test_point[1] = convert_to<Rational>(round_to_digits(d_test_point[1],9));
          rounded_test_point[2] = convert_to<Rational>(round_to_digits(d_test_point[2],9));
        
          rounded_test_point[0] = convert_to<Rational>(d_test_point[0]);
          rounded_test_point[1] =  convert_to<Rational>(d_test_point[1]);
          rounded_test_point[2] =  convert_to<Rational>(d_test_point[2]);

          //current_point = rounded_test_point;
         
        */
        current_point = test_point;

        cout << convert_to<double>(current_point) << endl;
        
        halfspace test_halfspace = minimum_halfspace_slices(slices_vertices, all_vertices, slices_triangulations, slices_edges,  TotalVolume, slices_volumes, slices_dimensions, ambient_dimension, centerpoint_slice_index, current_point, i);
        current_direction = - test_halfspace.direction;

         //Rescale the direction in the current slice to be unit norm

        Rational scale_denum_it = 1;

        for(int j = 0; j < ambient_dimension - centerpoint_slice_dimension; j++)
        {
            scale_denum_it -= current_direction[j]*current_direction[j];
        }

        Assert(scale_denum_it >= 0);

        if(scale_denum_it > 0)
        {
            for(int k = 0; k < current_direction.size(); k++)
            {
                //Maybe its better to leave the direction in double?
                if(k >= ambient_dimension - centerpoint_slice_dimension)
                {
                    double coordinate = convert_to<double>(current_direction[k]);
                    double new_coordinate = coordinate / std::sqrt(convert_to<double>(scale_denum_it));
                    current_direction[k] = convert_to<Rational>(new_coordinate);
                }
                else
                {
                    current_direction[k] = 0;
                }
            
                
            }
        }
        
       
    
   

        log_point *intermediate_point = PushStruct(Arena, log_point);
        intermediate_point->point.X = convert_to<real64>(current_point[0]);
        intermediate_point->point.Y = convert_to<real64>(current_point[1]);
        intermediate_point->point.Z = convert_to<real64>(current_point[2]);
        intermediate_point->direction.X = convert_to<real64>(current_direction[0]);
        intermediate_point->direction.Y = convert_to<real64>(current_direction[1]);
        intermediate_point->direction.Z = convert_to<real64>(current_direction[2]);
        intermediate_point->value = test_halfspace.value;
        intermediate_point->total_infeasibilities = test_halfspace.infeasibilities;

    }

    log_point best_found = *((log_point*) Arena->Base);

    for(int i = 0; i < (Arena->Used/sizeof(log_point)); i++)
    {
        log_point* step = (log_point *) (Arena->Base + i*sizeof(log_point));
        cout << "Iteration " << i << ": Point: " << step->point.X << " " << step->point.Y
             << " " << step->point.Z << "  -----> (next) direction: " << step->direction.X << " " << step->direction.Y << " " << step->direction.Z << " (value: " << step->value << ") " <<  " XXX: " << step->total_infeasibilities <<endl;

        if(step->value > best_found.value)
        {
            best_found = *step;
        }
     
    }

    cout << "(infeasabilities " << unfeasable_count << ")" << endl;

    cout << "Best " << best_found.value << endl;
    cout << "F( ) " << best_found.value / convert_to<double>(TotalVolume) << endl;

    return best_found;

}


void optimize_slices(Array<Matrix<Rational>> slices_vertices, Array<Vector<Rational>> points, memory_arena *Arena)
{
    int64 num_slices = slices_vertices.size();
    log_point slices_centerpoints[num_slices];
    
    for(int i = 0; i < num_slices; i++)
    {
        slices_centerpoints[i] = sub_grad_descend_slices(slices_vertices, points[i], i, Arena);
        //Reset Memory Arena:
        Arena->Used = 0;

    }

    log_point best_found = slices_centerpoints[0];
    for(int i = 0; i < num_slices; i++)
    {
        cout << "Slice " << i << ": Centerpoint: " <<  slices_centerpoints[i].point.X << " " << slices_centerpoints[i].point.Y << " " << slices_centerpoints[i].point.Z << "  -----> (next) direction: " << slices_centerpoints[i].direction.X << " " << slices_centerpoints[i].direction.Y << " " << slices_centerpoints[i].direction.Z << " (value: " << slices_centerpoints[i].value << ") "<< endl;

        if(slices_centerpoints[i].value > best_found.value)
        {
            best_found = slices_centerpoints[i];
        }

    }
    
    cout << endl;
    cout << "CENTERPOINT:  ( " << best_found.point.X << ", " << best_found.point.Y << ", " << best_found.point.Z << ")    VALUE: " << best_found.value <<  endl;


}
