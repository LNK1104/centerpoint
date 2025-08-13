halfspace minimum_halfspace_v2(BigObject* p, Matrix<Rational> vertices, Vector<Rational> centerpoint)
{


    BigObject simplicalcomplex;

    p->give("TRIANGULATION") >> simplicalcomplex;

    Array<Set<Int>> triangulation = simplicalcomplex.give("FACETS");

    BigObject simplicalcomplex_boundary;

    simplicalcomplex.give("BOUNDARY") >> simplicalcomplex_boundary;

     Array<Set<Int>> boundary_triangulation = simplicalcomplex_boundary.give("FACETS");


     int64 dimension = p->call_method("DIM");

     //Graph<> graph = p.give("GRAPH");
     Graph<> graph_adjacency = p->give("GRAPH.ADJACENCY");
     //IncidenceMatrix<> adjacency = p.give("GRAPH.EDGES");
     IncidenceMatrix<> adjacency = adjacency_matrix(graph_adjacency);
     BigObject graph("Graph<Undirected>");
     graph = p->give("GRAPH");
     //Graph<> BG=graph.give("ADJACENCY");
     Array<Set<Int>> edges = graph.call_method("EDGES");

      

    cout << boundary_triangulation << endl;

    Matrix<Rational> normals(vertices.rows(), vertices.cols());  // Create an empty matrix with correct size
    for (int i = 0; i < vertices.rows(); i++)
    {
        normals[i] = vertices.row(i) - centerpoint;  // Assign the same vector to each row
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
   
    Matrix<Rational> minimal_directions = zero_matrix<Rational>(signatures.rows(), dimension);
    
    
    for(int i = 0; i < signatures.rows(); i++)
    {
        cout << "-------------------------------------------------------------------------------------" << endl;
        cout << "Computing Halfspace Depth for Vertices: " << endl;

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

        /*
        Rational tentimes = 100000000000000;
        cone_h_representation /= tentimes;
        */

        Matrix<double> test_h_rep = convert_to<double>(cone_h_representation);
        cout << "++++++++++++++++++++++++++++++++" << endl;
        cout << test_h_rep << endl;
        cout << "++++++++++++++++++++++++++++++++" << endl;
        
       

        
        Matrix<Rational> sub_vertices;
        Matrix<Rational> conter_vertices;
        
       
        
        IncidenceMatrix<> total_adjacency = adjacency;
        IncidenceMatrix<> sub_total_adjacency = adjacency;

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
                sub_vertices /= vertices.row(j);
                sub_points.collect(j);
                sub_points_array.insert({sub_index, j});
                sub_index++;
            }
            else
            {
                conter_vertices /= vertices.row(j);
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

        Matrix<Int> edge_map = -1*ones_matrix<Int>(vertices.rows(), vertices.rows());

        
        //TODO: THAT doesn't work! Need the edges of the triangulation!!
     
        for(int j = 0; j < edges.size(); j++)
        {
            
           
            if((signatures(i,edges[j].front()) && !signatures(i,edges[j].back())))
            {
                
                sub_boundary_vertices /= vertices.row(edges[j].front());
                sub_boundary_points.collect(edges[j].front());
                sub_boundary_array.push_back(edges[j].front());

                conter_boundary_vertices /= vertices.row(edges[j].back());
                conter_boundary_points.collect(edges[j].back());
                conter_boundary_array.push_back(edges[j].back());

                edge_map(edges[j].front(), edges[j].back()) = number_lambda;
                number_lambda++;

                
                
            }
            else if((signatures(i,edges[j].back()) && !signatures(i,edges[j].front())))
            {
                sub_boundary_vertices /= vertices.row(edges[j].back());
                sub_boundary_points.collect(edges[j].back());
                sub_boundary_array.push_back(edges[j].back());
              
                
                conter_boundary_vertices /= vertices.row(edges[j].front());
                conter_boundary_points.collect(edges[j].front());
                conter_boundary_array.push_back(edges[j].front());

                edge_map(edges[j].back(), edges[j].front()) = number_lambda;
                number_lambda++;
                
            }
        }

        
       
        




        Matrix<Rational> qobj_coeff = zero_matrix<Rational>(number_lambda, number_lambda);
        Vector<Rational> lobj_coeff = zero_vector<Rational>(number_lambda);
        Rational FixVolume = 0;
        Rational ConterVolume = 0;
    
        


        for(int j = 0; j < boundary_triangulation.size(); j++)
        {

            std::vector<lambda_u_constraint> constraints;
            Vector<Int> triangle = Vector<Int>(boundary_triangulation[j]);


            for(int k = 0; k < triangle.size(); k++)
            {
                if(sub_boundary_points.contains(triangle[k]))
                {
                    for(int l = 0; l < triangle.size(); l++)
                    {
                        if(conter_boundary_points.contains(triangle[l]))
                        {
                            lambda_u_constraint cutting_edge = {};
                            cutting_edge.start_vec_index = triangle[k];
                            cutting_edge.end_vec_index = triangle[l];
                            constraints.push_back(cutting_edge);
                           
                        }

                    }
                       
                }
                
            }

            switch (dimension)
            {
                case 3:
                    // We Alwas have 0 or 2 contraits
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

                        


                            Vector<Rational> vec0 = vertices.row(constraints[0].start_vec_index);
                            Vector<Rational> vec1 = vertices.row(constraints[0].end_vec_index);
                            Vector<Rational> vec2 = vertices.row(constraints[1].end_vec_index);


                            //Computing Volume via determinant
                            Rational lambda_1_lambda_2 = cross_product3(vec1-vec0, vec2-vec0)*(centerpoint-vec0);

                        
                            qobj_coeff(edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index),
                                       edge_map(constraints[1].start_vec_index, constraints[1].end_vec_index)) =
                                Rational(1,6)*AbsoluteValue(lambda_1_lambda_2);
                               
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
                                Vector<Rational> vec0 = vertices.row(constraints[0].start_vec_index);
                            Vector<Rational> vec1 = vertices.row(constraints[1].start_vec_index);
                            Vector<Rational> vec2 = vertices.row(constraints[0].end_vec_index);


                            //Computing Volume via determinant
                            Rational lambda_1_lambda_2 = cross_product3(vec2-vec0, vec2-vec1)*(centerpoint-vec0);
                            Rational lambda_1 =cross_product3(vec2-vec0, vec1-vec0)*(centerpoint-vec0);

                            //Triangle 2: 
                            //Computing Volume via determinant
                       
                            Rational lambda_2 =cross_product3(vec2-vec1, vec1-vec0)*(centerpoint-vec0);
                            Rational constant = cross_product3(vec1-vec0, vec1-vec0)*(centerpoint-vec0);
                            Assert(convert_to<double>(constant) == 0);

                            if((cross_product3(vec2-vec0, vec1-vec0)*(centerpoint - vec0)) >= 0)
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
                        Vector<Rational> vec0 = vertices.row(triangle[0]);
                        Vector<Rational> vec1 = vertices.row(triangle[1]);
                        Vector<Rational> vec2 = vertices.row(triangle[2]);

                        Rational Volume = cross_product3(vec1-vec0, vec2-vec0)*(centerpoint-vec0);

                        if(sub_points.contains(triangle[0]))
                        {
                            FixVolume += Rational(1,6)*AbsoluteValue(Volume);
                        }
                        else
                        {
                            ConterVolume += Rational(1,6)*AbsoluteValue(Volume);
                        }

                    }
            }

        }

        //Build optimization problem

              
        bool32 clamp_precision = false;
        double    sol[number_lambda + dimension];
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
        int       qrow[dimension];
        int       qcol[dimension];
        double    qval[dimension];
        int       lind[dimension];
        double    lval[dimension];
        int       optimstatus;
       
        double    lambda_lb[number_lambda];
        double    lambda_ub[number_lambda];
        double    u_lb[dimension];
        double    u_ub[dimension];
        int       u2_qrow[dimension];
        int       u2_qcol[dimension];
        double    u2_qval[dimension];

        //NOTE: The bounds on the lambda are very important otherwise the computation takes forever... 
              for(int32 w = 0; w < number_lambda; w++)
          {
              lambda_lb[w] = 0;
              lambda_ub[w] = 1;
          }

           
        for(int32 w = 0; w < dimension; w++)
        {
            u_lb[w] = -1;
            u_ub[w] = 1;
        }
           
            
         

        for(int32 w = 0; w < lobj_coeff.size(); w++)
        {
            obj_lin[w] =  convert_to<double>(lobj_coeff[w]);
        }

            
            
                                      

        /* Create environment */

        error = GRBloadenv(&env, NULL);
        if (error) goto QUIT1;

#if OPT_NO_OUTPUT            
        // Stop gurobi from printing to the console
        error = GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0);
        if (error) goto QUIT1;
#endif

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
            error = GRBaddvars(model, dimension, 0, NULL, NULL, NULL, NULL, u_lb, u_ub, NULL, NULL);
            if (error) goto QUIT1;


            /* Add constraints */

            // Hyperplane constraint for each lambda
            //
            for(int32 k = 0; k < edge_map.rows(); k++)
            {
                for(int32 l = 0; l < edge_map.cols(); l++)
                {
                    Vector<Rational> vec0 = vertices.row(k);
                    Vector<Rational> vec1 = vertices.row(l);
                    Vector<Rational> num_vec = centerpoint - vec0;
                    Vector<Rational> denum_vec = vec1 - vec0;
                    
                    if(edge_map(k,l) != -1)
                    {
                        for(int32 v = 0; v < dimension; v++)
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
                
                        error = GRBaddqconstr(model, dimension, lind, lval, dimension, qrow, qcol, qval, GRB_EQUAL, 0, NULL);
                        if (error) goto QUIT1;
                    }
                }

            }
            
            
               

            // Cone constraint for u
            //
            //TODO: SCALE DOWN THE CONE VECTORS!!
            for(int32 w = 0; w < cone_h_representation.rows(); w++)
            {

                for(int32 v = 0; v < dimension; v++)
                {
                    //Note number_lambda + v is the u index
                    if(clamp_precision == true)
                    {
                        lind[v] = number_lambda+v; lval[v] = decimal_truncate(convert_to<double>(cone_h_representation(w,v)),3);
                        /*
                          cout << "lind " << v << " " << lval[v] << " " << convert_to<double>(test_cone(w,v)) << endl;
                          real64 truncate_test = truncate_to_digits_v2(convert_to<double>(test_cone(w,v)),3);
                          cout << truncate_test << endl;
                        */
                    }
                    else
                    {
                        lind[v] = number_lambda+v; lval[v] = convert_to<double>(cone_h_representation(w,v));
                    }
                
                }
                

                error = GRBaddconstr(model, dimension, lind, lval, GRB_GREATER_EQUAL, 0.0, NULL);
                if (error) goto QUIT1;

            }



            
            //Norm constraint for u, or just u > 0
            //u2_qrow[number_lambda] = 1.0;
            for(int32 w = 0; w < dimension; w++)
            {
                //Note number_lambda + w is the u index
                u2_qrow[w] = number_lambda+w; u2_qcol[w] = number_lambda+w; u2_qval[w]= 1.0;
            }
            
            error = GRBaddqconstr(model, 0, NULL, NULL, dimension, u2_qrow, u2_qcol, u2_qval, GRB_EQUAL, 1.0, "norm");
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
            
            
            
            
   
       
            

            /* Optimize model */

            error = GRBoptimize(model);
            if (error) goto QUIT1;

           

            /* Capture solution information */

            error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
            if (error) goto QUIT1;

            error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
            if (error) goto QUIT1;

            error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, number_lambda+dimension, sol);
            if (error) goto QUIT1;

            printf("\nOptimization complete\n");
            if (optimstatus == GRB_OPTIMAL) {
                printf("Optimal objective: %.4e\n", objval);

                for(int32 w = 0; w < number_lambda; w++)
                {
                    cout << "l" << w << "=" << sol[w]  << " ";
                   
                }

                //NOTE: WE cast the direction back to Rational; Do I have to worry about precision??
                
                for(int32 w = 0; w < dimension; w++)
                {                   
                    minimal_directions(i,w) = convert_to<Rational>(sol[number_lambda + w]);
                    cout << "u" << w << "=" << sol[number_lambda + w]  << " ";
                }

                
                
                cout << std::defaultfloat << std::setprecision(6);
                cout << endl;


                /* Free model */

                GRBwrite(model, "Test_positivity_qcp.lp");
           

                GRBfreemodel(model);

                /* Free environment */

                GRBfreeenv(env);

                //diffrent_contraint = true;
                //goto REPEAT;
               
                
            } else if (optimstatus == GRB_INF_OR_UNBD) {
                printf("Model is infeasible or unbounded\n");
            } else {
                printf("Optimization was stopped early\n");
            }

            
        
QUIT1:

            /* Error reporting */

            if (error) {
                printf("ERROR: %s\n", GRBgeterrormsg(env));
                

                error = GRBwrite(model, "error_qcp.lp");
                //Assert(false);

                if (error)
                {
                    cout << "Problem when writing error model" << endl;
                };

                if(clamp_precision == false)
                {

                    GRBfreemodel(model);
                    GRBfreeenv(env);
                
                    clamp_precision = true;
                    goto REPEAT;
                }
                else
                {
                    GRBwrite(model, "new_error_qcp.lp");
                    cout << "clambing didn't work" << endl;
                    exit(1);
                }
            }

        

            /* Free model */

            GRBfreemodel(model);

            /* Free environment */

            GRBfreeenv(env);

           

            Matrix<double> real_vertices = convert_to<double>(vertices);

            

           


            
            for(int j = 0; j < boundary_triangulation.size(); j++)
            {
                Vector<Int> triangle = Vector<Int>(boundary_triangulation[j]);
               
                std::vector<lambda_u_constraint> constraints;
                bool32 cut_triangle = false;

                for(int k = 0; k < triangle.size(); k++)
                {
                    if(sub_boundary_points.contains(triangle[k]))
                    {
                        for(int l = 0; l < triangle.size(); l++)
                        {
                            if(conter_boundary_points.contains(triangle[l]))
                            {
                                cut_triangle = true;
                                lambda_u_constraint cutting_edge = {};
                                cutting_edge.start_vec_index = triangle[k];
                                cutting_edge.end_vec_index = triangle[l];
                                constraints.push_back(cutting_edge);

                            }

                        }
                       
                    }

                   
                
                }

                if(cut_triangle == true)
                {
                    switch (dimension)
                    {
                        case 3:
                            if(constraints.size() > 0)
                            {
                                Assert(constraints.size() == 2);
                                cout << "Triangle: " << triangle << endl;
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

                        


                               

                                    double lambda1 = sol[edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index)];
                                    double lambda2 = sol[edge_map(constraints[1].start_vec_index, constraints[1].end_vec_index)];

                                    Vector<double> intersection_vec_1 = real_vertices.row(constraints[0].start_vec_index) + lambda1*(real_vertices.row(constraints[0].end_vec_index) - real_vertices.row(constraints[0].start_vec_index));
                                    Vector<double> intersection_vec_2 = real_vertices.row(constraints[1].start_vec_index) + lambda2*(real_vertices.row(constraints[1].end_vec_index) - real_vertices.row(constraints[1].start_vec_index));

                               
                                    Vector<double> vec1 = real_vertices.row(constraints[0].end_vec_index);
                                    Vector<double> vec2 = real_vertices.row(constraints[1].end_vec_index);

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
                                    real_conter_vertices /= convert_to<double>(centerpoint);
                                          
                                    
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

                                    double lambda1 = sol[edge_map(constraints[0].start_vec_index, constraints[0].end_vec_index)];
                                    double lambda2 = sol[edge_map(constraints[1].start_vec_index, constraints[1].end_vec_index)];

                                    Vector<double> intersection_vec_1 = real_vertices.row(constraints[0].start_vec_index) + lambda1*(real_vertices.row(constraints[0].end_vec_index) - real_vertices.row(constraints[0].start_vec_index));
                                    Vector<double> intersection_vec_2 = real_vertices.row(constraints[1].start_vec_index) + lambda2*(real_vertices.row(constraints[1].end_vec_index) - real_vertices.row(constraints[1].start_vec_index));

                               
                                    Vector<double> vec1 = real_vertices.row(constraints[0].end_vec_index);
                                   
                                    Matrix<double> real_conter_vertices = vector2row(intersection_vec_1);
                                    if((lambda1 != 1) || (lambda2 != 1))
                                    {
                                        real_conter_vertices /= intersection_vec_2;
                                        real_conter_vertices /= vec1;
                                    }
                                   
                                   
                                    real_conter_vertices /= convert_to<double>(centerpoint);
                                  
                                    
                                    Matrix<Rational> rational_conter_vertices = convert_to<Rational>(real_conter_vertices);

                                    BigObject conter_triangle("Polytope<Rational>");

                                    conter_triangle.take("VERTICES") << (ones_vector<Rational>() | rational_conter_vertices);
             
                                    Rational ConterVolume_triangle = conter_triangle.give("VOLUME");

                                    ConterVolume += ConterVolume_triangle;

                      


                        
                            
                                }
                    
                            }

                    }
                }
            }


           

            
            real64 sub_volume = objval + convert_to<double>(FixVolume);

            cout << "number_lambda " << number_lambda << endl;
            cout << "conter" << ConterVolume << endl;
           
            minimal_volume_signature_indices.push_back(i);
            minimal_volumes.push_back(sub_volume);
           
            
            real64 total_volume = convert_to<double>(ConterVolume) + objval + convert_to<double>(FixVolume);
            
            cout << "Volume: " << ConterVolume << "(conter) + " << objval << "(objval) + " << FixVolume << "(fix-sub) = " << total_volume << endl;


       
            if((total_volume > 167) || (total_volume < 166))
            {
                cout << "THE COMPUTED VOLUMES ARE NOT CORRECT !!!!";
                Assert(1 == 0);
                   
            }
        
       

        

       
       




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

    cout << "Minimal Direction: " << convert_to<double>(Halfspace.direction) << endl;

    return Halfspace;

    
    
}
