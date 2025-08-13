#include <polymake/Main.h>
#include <polymake/Matrix.h>
#include <polymake/SparseMatrix.h>
#include <polymake/Rational.h>
#include <polymake/IncidenceMatrix.h>
#include <polymake/Graph.h>
#include <polymake/client.h>
#include <polymake/Array.h>
#include <polymake/list>
#include <polymake/AccurateFloat.h>

#include <polymake/linalg.h>


#define OPT_DEBUG 0
#define OPT_NO_OUTPUT 0
#define SOL_GUROBI 1
#define SOL_BARON 0


#include "gurobi_c.h"


#include <fstream>
#include <sstream> // for iss in read_baron_results
#include <regex>

//for testing:
#include <random>



typedef int8_t int8;
typedef int16_t int16;
typedef int32_t int32;
typedef int64_t int64;

typedef uint8_t uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;

typedef int32 bool32;

typedef float real32;
typedef double real64;

/*
using namespace polymake;
int main(int argc, const char* argv[]) {
  try {
    const int dim = 4;
    Main pm;
    pm.set_application("polytope");
    BigObject p("Polytope<Rational>");
    
    p.take("VERTICES") << (ones_vector<Rational>() | 
       3*unit_matrix<Rational>(dim));
    const Matrix<Rational> f = p.give("FACETS");
    const Vector<Integer> h = p.give("H_STAR_VECTOR");
    cout << "facets" << endl << f << endl << "h* " << h << endl;
  } catch (const std::exception& ex) {
    std::cerr << "ERROR: " << ex.what() << endl; return 1;
  }
  return 0; 
}
*/

//TODO: make "get next conter element" a fuction
//TODO: make "add lambda" a function



#define Assert(Expression) if(!(Expression)){*(int *)0 = 0;}


#include "centerpoint_utilities.cpp"


using namespace polymake;

inline int64
RoundReal64ToInt64(real64 Real64)
{
    int64 Result = (int64)round(Real64);
    return Result;

}



inline real64 round_to_digits(real64 value, int64 digits)
{

    real64 Result = RoundReal64ToInt64(value * (pow(10, digits))) / (real64)(pow(10, digits));
    return  Result;
}

inline real64 decimal_truncate(real64 value, int64 digits)
{

    // cout <<  trunc(value * 10 * digits)  << endl;
    real64 Result = value / (real64)(pow(10, digits));
    return  Result;
}



struct lambda_u_constraint
{
    Vector<Rational> num_vec;
    Vector<Rational> denum_vec;

    //For debugging:
    int64 start_vec_index;
    int64 end_vec_index;
};

struct halfspace
{
    Vector<Rational> direction;
    real64 value;
    int infeasibilities;
};

/*
struct Pair
{
    int64 first;
    int64 second;
};
*/

struct Vec3
{
    real64 X;
    real64 Y;
    real64 Z;

};

struct log_point
{
    Vec3 point;
    Vec3 direction;
    real64 value;
    int total_infeasibilities;
};

struct interior_point
{
    /*
    Vector<Rational> point;
    Vector<Rational> direction;
    */
    real64 point_x;
    real64 point_y;
    real64 point_z;
    real64 direction_x;
    real64 direction_y;
    real64 direction_z;

    real64 value;
    real64 volume_difference;
    real64 step_size;
    bool32 found;
    real64 line_values[100];
    Vec3 line_points[100];
    Vec3 line_directions[100];
    int32 path;
};



inline Rational
AbsoluteValue(Rational value)
{
    //real32 Result = fabsf(Real32);
    Rational Result;
    if(value >= 0)
    {
        Result = value;
    }
    else
    {
        Result = -value; 
    }
    return(Result);
}



bool32 InFacet(int64 i, int64 j,const pm::IncidenceMatrix<> * facets)
{
    bool32 value = false;

    for(int l = 0; l < facets->rows(); l++)
    {
        if((*facets)(l, i) && (*facets)(l,j))
        {
            value = true;
        }
    }

    return value;
}




                
inline int64 getnextElement(int64 current, int64 next, std::vector<long> boundary_array,  IncidenceMatrix<> adjacency)
{

    long privious = current;
    current = next;
    int64 next_elements = 0;
                

               
    //Get next element
    if(boundary_array.size() == 2)
    {
        next = privious;
    }
    else{
        for(int l = 0; l < boundary_array.size(); l++)
        {
            if((adjacency(current, boundary_array[l])) && (boundary_array[l] != privious))
            {
                //NOTE: This path should bet taken only by one unique element!! No break statement!
                next = boundary_array[l];
                next_elements++;
            }    
        }

        Assert(next_elements <= 1);

    }

    return next;

}

inline real64 correctVolume(std::vector<long> sub_boundary_array,  IncidenceMatrix<> sub_adjacency, std::unordered_map<int, int>  sub_points_array, IncidenceMatrix<> sub_total_adjacency, IncidenceMatrix<> total_adjacency, Matrix<Rational> vertices, Vector<Rational> centerpoint)
{
    if(sub_boundary_array.size() > 3)
    {
        int64 current = sub_boundary_array[0];
        int64 next = current;

        for(int64 i = 0; i < sub_boundary_array.size(); i++)
        {
            if(total_adjacency(current, sub_boundary_array[i]))
            {
                next = sub_boundary_array[i];
            }
        }
        
        
        //Walk along the sub boundary
        for(int64 i = 0; i < sub_boundary_array.size(); i++)
        {
            for(int64 j = 0; j < sub_boundary_array.size(); j++)
            {
                if(sub_total_adjacency(current,sub_boundary_array[j]) && !total_adjacency(current,sub_boundary_array[j]) )
                {
                    //SUB CYCLE IN sub boundary
                    
                    cout << "CYCLE DETECTED!" << "( " << current << " - " << sub_boundary_array[j] << " )" <<  endl;

                    //check if the next element is connected with the centerpoint
                    bool32 connected_with_centerpoint = false;
                    for(int k = 0; k < sub_adjacency.cols()-1; k++)
                    {
                        if((sub_points_array[k] == next) && sub_adjacency(sub_adjacency.cols()-1, k))
                        {
                            connected_with_centerpoint = true;
                        }
                    }

                    if(connected_with_centerpoint == false)
                    {

                        Matrix<Rational> cycle_vertices;
                    

                        while(current != sub_boundary_array[j])
                        {
                            cycle_vertices /= vertices.row(current);
                            cout << current << " - ";  
               
                            int64 next_current = next;
                            next = getnextElement(current, next, sub_boundary_array, total_adjacency);
                            current= next_current;
                        }

                        cycle_vertices /= vertices.row(sub_boundary_array[j]);
                        cycle_vertices /= centerpoint;

                        BigObject subp("Polytope<Rational>");

                        subp.take("VERTICES") << (ones_vector<Rational>() | cycle_vertices);
        
                        real64 SurplusVolume = subp.give("VOLUME");

                        return SurplusVolume;
                    }
                }
            }



            //Get next element
            int64 next_current = next;
            next = getnextElement(current, next, sub_boundary_array, total_adjacency);
            current= next_current;
        }
                   
        
    }
    

    return 0;
}


Vector<Rational> cross_product3(Vector<Rational> vec1, Vector<Rational> vec2)
{
    Vector<Rational> solution = zero_vector<Rational>(3);
    

    solution[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
    solution[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
    solution[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
    

    return solution;
}

bool32 is_feasible(Matrix<Rational>* Inequalities, Vector<Rational> point)
{
    bool feasible = true;

    Vector<Rational> rhs = (*Inequalities)*(ones_vector<Rational>(1) | point);

    for(int64 i = 0; i < rhs.dim(); i++)
    {
        if(rhs[i] < 0)
        {
            feasible = false;
        }
    }

    return feasible;
}


#include "centerpoint_v2.cpp"
#include "centerpoint_slices.cpp"


halfspace minimum_halfspace(BigObject* p, Matrix<Rational> vertices, Vector<Rational> centerpoint)
{

    const IncidenceMatrix<> facets = p->give("VERTICES_IN_FACETS");

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
   


    //Graph<> graph = p.give("GRAPH");
    Graph<> graph_adjacency = p->give("GRAPH.ADJACENCY");
    //IncidenceMatrix<> adjacency = p.give("GRAPH.EDGES");
    IncidenceMatrix<> adjacency = adjacency_matrix(graph_adjacency);
    BigObject graph("Graph<Undirected>");
    graph = p->give("GRAPH");
    //Graph<> BG=graph.give("ADJACENCY");
    Array<Set<Int>> edges = graph.call_method("EDGES");
    cout << edges << endl;

    // graph = p.give("GRAPH");
    
    //Array<Set<Int>> edges2 = graph.call_method("EDGES");

    //cout << edges2 << endl;

    
    

    
   
        //const Vector<Integer> h = p.give("H_STAR_VECTOR");
    //cout << "facets" << endl << facets << endl << "h* " << adjacency << endl;
    cout << all_rays << endl;

    //Loop over all signaurtes
    cout << signatures << endl;

   
    std::vector<int64> minimal_volume_signature_indices;
    std::vector<real64> minimal_volumes;
    int32 dimension = 3;
    Matrix<Rational> minimal_directions = zero_matrix<Rational>(signatures.rows(), dimension);
    
    
    for(int i = 0; i < signatures.rows(); i++)
    {
        cout << "-------------------------------------------------------------------------------------" << endl;
        cout << "Computing Halfspace Depth for Vertices: ";
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

        //TODO: Improve Hash map!!! --- or use vector instead??


        
        //TEST: Compute edges that paricipate in cut:

       
        //Not working:
        //Set<Integer> chamber = HA.call_method("signature_to_chamber", signatures.row(i));

        //For Signatures: The i-th entry is the signature of the i-th maximal cone of the CHAMBER_DECOMPOSITION
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



        /*_________sub_polytope___________*/


        sub_vertices /= centerpoint;
        BigObject subp("Polytope<Rational>");

        subp.take("VERTICES") << (ones_vector<Rational>() | sub_vertices);
        //Matrix<Rational> sub_V = subp.give("VERTICES");
        //cout << "Subpolytope Vertices: " << sub_V << endl;
        //ASSERT: NO POINTS should be removed from the convex hull (ecept the centerpoint maybe...)
        //TODO: Seperatly check the edge cases (where lambda is 0)
        Graph<> sub_graph_adjacency = subp.give("GRAPH.ADJACENCY");
        IncidenceMatrix<> sub_adjacency = adjacency_matrix(sub_graph_adjacency);

        // TODO: WE CHEAT INSTEAD!!!!!

        //NOTE: Last vertex is centerpoint
        for(int j= 0; j < sub_adjacency.rows()-1; j++)
        {
            //NOTE: Last vertex is centerpoint
            for(int k= 0; k < sub_adjacency.cols()-1; k++)
            {
                if(InFacet(sub_points_array[j], sub_points_array[k], &facets))
                {
                    total_adjacency(sub_points_array[j], sub_points_array[k]) = sub_adjacency(j,k);
                }
            }
        }

        //NOTE: Last vertex is centerpoint
        for(int j= 0; j < sub_adjacency.rows()-1; j++)
        {
            //NOTE: Last vertex is centerpoint
            for(int k= 0; k < sub_adjacency.cols()-1; k++)
            {
               
                sub_total_adjacency(sub_points_array[j], sub_points_array[k]) = sub_adjacency(j,k);
                
            }
        }

        

        //NOTE: For lower dimensional Polytopes the Volume is 0!
        real64 FixVolume = subp.give("VOLUME");
        int64 dim = subp.give("CONE_DIM");

         //NOTE: The cone Dimension is "One more than the dimension of the affine hull of the polyhedron"
        dim = dim - 1;
        
        cout << "FixVolume: " << FixVolume << endl;
        cout << "Dimension of subp: " << dim << endl;

        

        /*_________conter_polytope___________*/
        

        //TODO Assert: Set<Int> conter_points = points - sub_points - centerpoint;



       

        //NOTE: It is necessary to add the centerpoint here, since otherwise the surceface pointing to the centerpoint gets connected via diagonals!!!

        
        conter_vertices /= centerpoint;

        //NOTE: Need to compute the adjacency matritx of the conter_polytope to ensure connectivity of "bourdrary" next to cutting volume

        // Matrix<Rational> conter_vertices = vertices.minor(conter_points, All);
        BigObject conterp("Polytope<Rational>");
        conterp.take("VERTICES") << (ones_vector<Rational>() | conter_vertices);

        Graph<> conter_graph_adjacency = conterp.give("GRAPH.ADJACENCY");
        IncidenceMatrix<> conter_adjacency = adjacency_matrix(conter_graph_adjacency);
       

        cout << total_adjacency << endl;
        
     
        
        //TODO WE Cheat instead
        
        for(int j= 0; j < conter_adjacency.rows()-1; j++)
        {
            for(int k= 0; k < conter_adjacency.cols()-1; k++)
            {
                if(InFacet(conter_points_array[j], conter_points_array[k], &facets))
                {
                    total_adjacency(conter_points_array[j], conter_points_array[k]) = conter_adjacency(j,k);
                }
            }
        }
        


        int32 number_lambda = 0;


        
        Set<Int> sub_boundary_points;
     
        Set<Int> conter_boundary_points;
     
        for(int j = 0; j < edges.size(); j++)
        {
            //NOTE: The edgelist alwas has two elements
            if((signatures(i,edges[j].front()) && !signatures(i,edges[j].back())))
            {
                number_lambda++;
                sub_boundary_points.collect(edges[j].front());
               
                conter_boundary_points.collect(edges[j].back());
                
                
            }
            else if((signatures(i,edges[j].back()) && !signatures(i,edges[j].front())))
            {
                number_lambda++;
                sub_boundary_points.collect(edges[j].back());
                
                conter_boundary_points.collect(edges[j].front());
                
            }
        }

        cout << "sub boundary: " << sub_boundary_points << endl;
        cout << "conter_boundary: " << conter_boundary_points << endl;

        cout << "Lambdas: " << number_lambda << endl;
        cout << total_adjacency << endl;
        cout << sub_total_adjacency << endl;

        




        Matrix<Rational> qobj_coeff = zero_matrix<Rational>(number_lambda, number_lambda);
        Vector<Rational> lobj_coeff = zero_vector<Rational>(number_lambda);
        Rational cobj_coeff = 0;
        std::vector<lambda_u_constraint> constraints;

       
        //TEST: check degree 2 connectivity of every point on the sub_boundary and conter_boundary

        if(sub_boundary_points.size() > 2)
        {
            for (long elem :sub_boundary_points)
            {
                int degree = 0;
                for (long elem2 :sub_boundary_points)
                {
                     if(total_adjacency(elem, elem2))
                    {
                        degree++;
                    }
                }

                if(degree != 2)
                {
                    cout << "ATTTENTION!!!! degree= " << degree << endl;
                }
            }

        }

        if(conter_boundary_points.size() > 2)
        {
            for (long elem :conter_boundary_points)
            {
                int degree = 0;
                for (long elem2 :conter_boundary_points)
                {
                     if(total_adjacency(elem, elem2))
                    {
                        degree++;
                    }
                }

                if(degree != 2)
                {
                    cout << "ATTTENTION!!!! degree= " << degree << endl;
                    
                }
            }

        }

       

        //Consturct sub_boudary_array/conter boundary array

        std::vector<long> sub_boundary_array = {};

        for(long elem :sub_boundary_points)
        {
            sub_boundary_array.push_back(elem);
        }

        std::vector<long> conter_boundary_array = {};

        for(long elem :conter_boundary_points)
        {
            conter_boundary_array.push_back(elem);
        }


        //FixVolume = correctVolume(int64 current, int64 next, std::vector<long>* sub_boundary_array,  IncidenceMatrix<>* total_adjacency)

        real64 differenceVolume = correctVolume(sub_boundary_array, sub_adjacency, sub_points_array, sub_total_adjacency, total_adjacency, vertices, centerpoint);

        cout << "Volume correction of Value: " << differenceVolume << endl;
        FixVolume -= differenceVolume; 


        //IDEA: For the 3 dimensional case: instead of hunting for the corresponding wedges for the cutting volume, walk the sub- and conter bundary... More complex in heigher dimensions


       
        long current = sub_boundary_array[0];
        long next = sub_boundary_array[0];
        long current_conter = 0;
        long next_conter = 0;

        int64 lambda_index = 0;
        int64 first_sub = current;
        int64 first_conter = 0;

        bool32 move_in_conter = false;
        bool32 reverse_dir = false;


       
        if(number_lambda == 5)
        {
            int testw = 0;
        }
        

        for(int k = 0; k < sub_boundary_array.size(); k++)
        {
            if(k == 0)
            {
                //Determines walk direction

                
                //Get first adjacent, if exists
                
                for(int l = 0; l < sub_boundary_array.size(); l++)
                {
                    if(total_adjacency(sub_boundary_array[0], sub_boundary_array[l]))
                    {
                        next = sub_boundary_array[l];
                        break;
                    }
                    
                        
                }

                
                int64 current_neighbour;

                //NOTE: Determine neighbour exactly as in the case down below
                //TODO: Make this one instead of two seperate... but we might have to recompute

                Set<Int> pre_adjacent_conter_to_current = {};
                std::vector<int64> pre_adjacent_conter_to_current_array = {};
                Set<Int> pre_adjacent_conter_to_next = {};
                std::vector<int64> pre_adjacent_conter_to_next_array = {};

                for(int l = 0; l < conter_boundary_array.size(); l++)
                {
                    if(total_adjacency(current, conter_boundary_array[l]))
                    {
                        pre_adjacent_conter_to_current.collect(conter_boundary_array[l]);
                        pre_adjacent_conter_to_current_array.push_back(conter_boundary_array[l]);
                    }
                }

                for(int l = 0; l < conter_boundary_array.size(); l++)
                {
                    if(total_adjacency(next, conter_boundary_array[l]))
                    {
                        pre_adjacent_conter_to_next.collect(conter_boundary_array[l]);
                        pre_adjacent_conter_to_next_array.push_back(conter_boundary_array[l]);
                    }
                }

                             
                Set<Int> pre_intersection_point = pre_adjacent_conter_to_current * pre_adjacent_conter_to_next;

                if(pre_intersection_point.size() >= 1)
                {
                    //NOTE: normal Traingle doesn't matter since they do not have a lambda 1!
                        
                    //degenerate triangle
                    current_neighbour = pre_intersection_point.front();
                }
                else
                {
                    //rectangle
                    Set<Int> connection_point_current = {};
                    Set<Int> connection_point_next = {};

                    for(int l = 0; l < pre_adjacent_conter_to_current_array.size(); l++)
                    {
                        for(int m = 0; m < pre_adjacent_conter_to_next_array.size(); m++)
                        {
                          
                            if(total_adjacency(pre_adjacent_conter_to_current_array[l], pre_adjacent_conter_to_next_array[m]))
                            {
                                connection_point_next.collect(pre_adjacent_conter_to_next_array[m]);
                                connection_point_current.collect(pre_adjacent_conter_to_current_array[l]);
                                break;
                            }
                        }
                    }


                    current_neighbour = connection_point_current.front();

                }

               

                Vector<Rational> vecx1 = vertices.row(next) - vertices.row(current);
                Vector<Rational> vecy2 = vertices.row(current_neighbour) -vertices.row(current);
                Vector<Rational> vecz = centerpoint - vertices.row(current);

                //NOTE: THE order vec2 and then vec1 is important.
                //This must coinside with ter order in the computation of lambda1 in the (degen triangle) and rectangle case!
                if((cross_product3(vecy2, vecx1)*vecz) < 0)
                {
                    //reverse the direction:

                    cout << "reverse direction! " << endl;

                    reverse_dir = true;

                    //NOTE: The case where the sub_boundary is only on point flows through the else case

                    if(sub_boundary_array.size() == 2)
                    {
                        int64 old_current = current;
                        current = next;
                        next = old_current;
                        first_sub = current;
                    }
                    else
                    {

                        for(int l = 0; l < sub_boundary_array.size(); l++)
                        {
                            if((total_adjacency(sub_boundary_array[0], sub_boundary_array[l])) && (sub_boundary_array[l] != next))
                            {
                                next = sub_boundary_array[l];
                                break;
                            }
                    
                        
                        }
                    }

                    
                }
                


                cout << "direction: " << current << ", " << next << endl;


                // ----------- end of determining walk direction ---------------




                

                //Determin all adjacent conter_points
                //Path if sub_polytope contains only one point

                //Sync conter_boudary direction with choosen sub_boundary direction

              
                Set<Int> adjacent_conter_to_current = {};
                std::vector<int64> adjacent_conter_to_current_array = {};
                Set<Int> adjacent_conter_to_next = {};
                std::vector<int64> adjacent_conter_to_next_array = {};
                

                for(int l = 0; l < conter_boundary_array.size(); l++)
                {
                    if(total_adjacency(current, conter_boundary_array[l]))
                    {
                        adjacent_conter_to_current.collect(conter_boundary_array[l]);
                        adjacent_conter_to_current_array.push_back(conter_boundary_array[l]);
                    }
                }

                for(int l = 0; l < conter_boundary_array.size(); l++)
                {
                    if(total_adjacency(next, conter_boundary_array[l]))
                    {
                        adjacent_conter_to_next.collect(conter_boundary_array[l]);
                        adjacent_conter_to_next_array.push_back(conter_boundary_array[l]);
                    }
                }

                


                if(sub_points.size() == 1)
                {
                    //triangle

                    /* case:
               current |----------| adjacent_conter_to_current_array[l]
                       | \        |
                       |    \     |
                       |      \   |
                       |         \|
                       |----------| adjacent_conter_to_current_array[m]

                    */

                    //Determine walk direction on the conterboundary
                    bool32 found = false;
                    for(int l = 0; l < adjacent_conter_to_current_array.size(); l++)
                    {
                        for(int m = 0; m < adjacent_conter_to_current_array.size(); m++)
                        {
                            if(total_adjacency(adjacent_conter_to_current_array[l], adjacent_conter_to_current_array[m]))
                            {
                                //NOTE: the two points should be in the same facet
                                Assert(InFacet(adjacent_conter_to_current_array[l], adjacent_conter_to_current_array[m], &facets));
                                current_conter = adjacent_conter_to_current_array[l];
                                next_conter = adjacent_conter_to_current_array[m];
                                first_conter = current_conter;
                                found = true;
                                break;
                                    
                            }
                        }
                        if(found)
                        {
                            break;
                        }
                    }


                    //Do the Triangles in chosen order
                    do
                    {
                        Vector<Rational> vec0 = vertices.row(current);
                        Vector<Rational> vec1 = vertices.row(current_conter);
                        Vector<Rational> vec2 = vertices.row(next_conter);


                        //Computing Volume via determinant
                        Rational lambda_1_lambda_2 = cross_product3(vec1-vec0, vec2-vec0)*(centerpoint-vec0);
                        qobj_coeff(lambda_index, (lambda_index +1)%number_lambda) = Rational(1,6)*AbsoluteValue(lambda_1_lambda_2);
                               

                        //NOTE:
                        //
                        //vec0 is in the subboundary and vec1 in the conterboundary
                        //lambda is defined via vec0 + lambda*(vec1-vec0) = lambda*vec1 + (1-lambda)vec0
                        //We are always pushing the lambda of the upper line!!
                        //
                        //
                        Vector<Rational> l1_nominator= centerpoint-vec0;
                        Vector<Rational> l1_denominator = vec1-vec0;

                        
                        constraints.push_back({l1_nominator, l1_denominator, current, current_conter});
                        lambda_index++;

                        

                        cout << "first conter (Triangle): " << current_conter << ", " << next_conter << endl;


                        
                        //Get next conter element

                        long next_next_conter = next_conter;
                        
                        if(conter_boundary_array.size() == 2)
                        {
                            if(next_conter == conter_boundary_array[0])
                            {
                                next_next_conter = conter_boundary_array[1];
                            }
                            else
                            {
                                next_next_conter = conter_boundary_array[0];
                            }
                        }
                        else
                        { 
                            int64 next_conter_elements = 0;
                            for(int m = 0; m < conter_boundary_array.size(); m++)
                            {
                                if((total_adjacency(next_conter, conter_boundary_array[m])) && (conter_boundary_array[m] != current_conter))
                                {
                                    //NOTE: This path should bet taken only by one unique element!! No break statement!
                                    next_next_conter = conter_boundary_array[m];
                                    next_conter_elements++;
                                }    
                            }

                            //NOTE: NO next elment if subboundary has two elements...
                     
                            Assert(next_conter_elements <= 1);
                        }

                        current_conter = next_conter;
                        next_conter = next_next_conter;
                        
                    }
                    while(current_conter != first_conter);
                    

                    
                }

                else
                {

                    //Check for triangle with next point
                                
                    Set<Int> intersection_point = adjacent_conter_to_current * adjacent_conter_to_next;

                    //NOTE: If we have more intersection_points (e.g. for a 3-dimensional simplex) we just use the fist one
                    if(intersection_point.size() >= 1)
                    {
                        /* case:
                   current |----------| intersection_index
                           |        / |
                           |      /   |
                           |   /      |
                           | /        |
                      next |----------|

                        */

                        long intersection_index = intersection_point.front();

                        //dimension3
                    
                        if(vertices.cols() == 3)
                        {

                            //Triangle 1: 
                                Vector<Rational> vec0 = vertices.row(current);
                            Vector<Rational> vec1 = vertices.row(next);
                            Vector<Rational> vec2 = vertices.row(intersection_index);


                            //Computing Volume via determinant
                            Rational lambda_1_lambda_2 = cross_product3(vec2-vec0, vec2-vec1)*(centerpoint-vec0);
                            Rational lambda_1 =cross_product3(vec2-vec0, vec1-vec0)*(centerpoint-vec0);

#if OPT_DEBUG
                        
                        
                            cout << "l1_l2: " << convert_to<double>(lambda_1_lambda_2)  << endl;
                            cout << "l1: " << convert_to<double>(lambda_1) << endl;

                            cout << convert_to<double>(qobj_coeff) << endl;
                            cout << convert_to<double>(lobj_coeff) << endl;
                          
                   
                        

                        

#endif
                         

                            //NOTE in lobj_coeff += necessary since we might push linear terms ahead
                          
                            if(lambda_1_lambda_2 + lambda_1 < 0)
                            {
                                qobj_coeff(lambda_index, (lambda_index +1)%number_lambda) = Rational(1,6)*(-lambda_1_lambda_2);
                                lobj_coeff[lambda_index] += Rational(1,6)*(-lambda_1);
                            
                            }
                            else
                            {
                                qobj_coeff(lambda_index, (lambda_index +1)%number_lambda) = Rational(1,6)*(lambda_1_lambda_2);
                                lobj_coeff[lambda_index] += Rational(1,6)*(lambda_1);
                            }
                          
                               
                                
                            Vector<Rational> l1_nominator= centerpoint-vec0;
                            Vector<Rational> l1_denominator = vec2-vec0;
                            constraints.push_back({l1_nominator, l1_denominator, current, intersection_index});
                        
                        
                            /*
                              Matrix<Rational> l1_l2_u_nominator = zero_matrix<Rational>(vertices.cols(), vertices.cols());
                              for(int n = 0; n < vertices.cols(); n++)
                              {
                              l1_l2_u_nominator[n] = lambda_1_lambda_2*(centerpoint - vec0)[n]*(centerpoint-vec1);
                              }

                              Matrix<Rational> l1_l2_u_denominator = zero_matrix<Rational>(vertices.cols(), vertices.cols());
                              for(int n = 0; n < vertices.cols(); n++)
                              {
                              l1_l2_u_denominator[n] = (vec0-vec2)[n]*(vec1-vec2);
                              }

                              Vector<Rational> l1_u_nominator = lambda_1*(centerpoint-vec0);
                              Vector<Rational> l1_u_denominator = (vec0-vec2);
                            */

                            //Triangle 2: 
                       

                            //Computing Volume via determinant
                       
                            Rational lambda_2 =cross_product3(vec1-vec0, vec2-vec1)*(centerpoint-vec0);
                            Rational constant = cross_product3(vec1-vec0, vec1-vec0)*(centerpoint-vec0);

#if OPT_DEBUG
                        
                        
                          cout << "l2: " << lambda_2  << endl;
                          cout << "l1: " << constant << endl;

                          cout << qobj_coeff << endl;
                          cout << lobj_coeff << endl;
                          
                   
                        

                        

#endif
                          
      
                            if(lambda_2 + constant < 0)
                            {
                           
                                lobj_coeff[(lambda_index +1)%number_lambda] += Rational(1,6)*(-lambda_2);
                                cobj_coeff += Rational(1,6)*(-constant);
                            
                            }
                            else
                            {
                                lobj_coeff[(lambda_index +1)%number_lambda] += Rational(1,6)*(lambda_2);
                                cobj_coeff += Rational(1,6)*(constant);
                            }

                          
                        
                            lambda_index++;
                                    
                        }
                    
                        first_conter = intersection_index;
                        cout << "first conter (degenerate Triangle): " << intersection_index << endl;

                        //___Determine the conter points to follow_____

                        if(conter_boundary_array.size() == 1)
                        {
                            current_conter = intersection_index;
                            next_conter = intersection_index;
                        }
                        else if(conter_boundary_array.size() == 2)
                        {
                            if(intersection_index == conter_boundary_array[0])
                            {
                                current_conter = intersection_index;
                                next_conter = conter_boundary_array[1];
                            }
                            else
                            {
                                current_conter = intersection_index;
                                next_conter = conter_boundary_array[0];
                            }
                        }
                        else
                        {

                            // Determine the privious conter_index to sync subbourdary direction with conterboudary direction

                            //follow subbourdrary backward until we are at a point adjacent to intersection index neighbors;

                            current_conter = intersection_index;
                            long intersection_index_neighbor1 = intersection_index;
                            long intersection_index_neighbor2 = intersection_index;

                            int64 traversal_neighbors = 0;

                            for (int m = 0; m < conter_boundary_array.size(); m++)
                            {
                                if(total_adjacency(intersection_index, conter_boundary_array[m]))
                                {
                                    if(traversal_neighbors == 0)
                                    {
                                        intersection_index_neighbor1 = conter_boundary_array[m];
                                        traversal_neighbors++;
                                    }
                                    else
                                    {
                                        intersection_index_neighbor2 = conter_boundary_array[m];
                                        traversal_neighbors++;
                                    }
                                
                                }
                            }

                       

                            //Compute adjacent vertices to neighors in subboundary
                            Set<Int> adjacent_to_neighbor1 = {};
                            Set<Int> adjacent_to_neighbor2 = {};


                            for (int m = 0; m < sub_boundary_array.size(); m++)
                            {
                                if(total_adjacency(intersection_index_neighbor1, sub_boundary_array[m]))
                                {
                                    adjacent_to_neighbor1.collect(sub_boundary_array[m]);
                                
                                }
                                if(total_adjacency(intersection_index_neighbor2, sub_boundary_array[m]))
                                {
                                    adjacent_to_neighbor2.collect(sub_boundary_array[m]);
                                }
                        
                            }

                            //Walk along subboundary in next direction to see, which neighbors adjacencys get visited first
                            long privious_lookahead = current;
                            long current_lookahead = next;
                            for(int m = 0; m < sub_boundary_array.size(); m++)
                            {
                                
                               
                                //NOTE: we need to consider this case: see ex {1,2,3,7}
                                if(adjacent_to_neighbor1.contains(current_lookahead) && adjacent_to_neighbor2.contains(current_lookahead))
                                {
                                    if(!total_adjacency(intersection_index_neighbor1, privious_lookahead))
                                    {
                                        current_conter = intersection_index;
                                        next_conter = intersection_index_neighbor1;
                                        break;
                                    }
                                    else
                                    {
                                        current_conter = intersection_index;
                                        next_conter = intersection_index_neighbor2;
                                        break;
                                    }

                                         
                                }
                                else{
                                        
                                    if(adjacent_to_neighbor1.contains(current_lookahead))
                                    {
                                        current_conter = intersection_index;
                                        next_conter = intersection_index_neighbor1;
                                        break;
                                    
                                    }
                                    else if(adjacent_to_neighbor2.contains(current_lookahead))
                                    {
                                        current_conter = intersection_index;
                                        next_conter = intersection_index_neighbor2;
                                        break;
                                        
                                    }
                                }

                                        

                                long next_lookahead = current_lookahead;
                                int64 next_elements = 0;
                

                                
                                //Get next element
                                if(dim == 2)
                                {
                                    next_lookahead = privious_lookahead;
                                }
                                else{
                                    for(int l = 0; l < sub_boundary_array.size(); l++)
                                    {
                                        if((total_adjacency(current_lookahead, sub_boundary_array[l])) && (sub_boundary_array[l] != privious_lookahead))
                                        {
                                            //NOTE: This path should bet taken only by one unique element!! No break statement!
                                            next_lookahead = sub_boundary_array[l];
                                            next_elements++;
                                        }    
                                    }

                                    Assert(next_elements <= 1);

                                }

                                privious_lookahead = current_lookahead;
                                current_lookahead = next_lookahead;

                           
                            }
                        
                        }
                    }

                    //NOTE: the case where intersection_point.size() is not 0 or 1 is only when subpolytope is 1 point --- this is handeled already in first case
                    else if(intersection_point.size() == 0)
                    {
                        //Rectangle: Triangulate
                    

                    
                        /* case:
                   current |----------| connection_index_current
                           |          |
                           |          |
                           |          |
                           |          |
                      next |----------| connection_index_next

                        */
                    

                        //Search for connetion in the two sets.
                   
                        Set<Int> connection_point_current = {};
                        Set<Int> connection_point_next = {};

                        for(int l = 0; l < adjacent_conter_to_current_array.size(); l++)
                        {
                            for(int m = 0; m < adjacent_conter_to_next_array.size(); m++)
                            {
                          
                                if(total_adjacency(adjacent_conter_to_current_array[l], adjacent_conter_to_next_array[m]))
                                {
                                    connection_point_next.collect(adjacent_conter_to_next_array[m]);
                                    connection_point_current.collect(adjacent_conter_to_current_array[l]);
                                    break;
                                }
                            }
                        }

                        //NOTE: The break in the upper if is only neccessary, since in the case of sub_point.size()                    
                                                              
                        long connection_index_next = connection_point_next.front();
                        long connection_index_current = connection_point_current.front();


                     
                        Vector<Rational> vec0 = vertices.row(current);
                        Vector<Rational> vec1 = vertices.row(next);
                        Vector<Rational> vec2 = vertices.row(connection_index_current);
                        Vector<Rational> vec3 = vertices.row(connection_index_next);
                   

                        //Triangle 1: 
                        //Computing Volume via determinant

                     
                        Rational lambda_1_lambda_2 = cross_product3(vec2-vec0, vec3-vec1)*(centerpoint-vec0);
                        Rational lambda_1 = cross_product3(vec2-vec0, vec1-vec0)*(centerpoint-vec0);
#if OPT_DEBUG
                        
                        
                          cout << "l1_l2: " << lambda_1_lambda_2  << endl;
                          cout << "l1: " << lambda_1 << endl;
                          
                          cout << qobj_coeff << endl;
                          cout << lobj_coeff << endl;
                          
                   
                        

                        

#endif
                          

                        if(lambda_1_lambda_2 + lambda_1 < 0)
                        {
                            qobj_coeff(lambda_index, (lambda_index +1)%number_lambda) = Rational(1,6)*(-lambda_1_lambda_2);
                            lobj_coeff[lambda_index] += Rational(1,6)*(-lambda_1);
                            
                        }
                        else
                        {
                            qobj_coeff(lambda_index,(lambda_index +1)%number_lambda) = Rational(1,6)*(lambda_1_lambda_2);
                            lobj_coeff[lambda_index] += Rational(1,6)*(lambda_1);
                        }
                          

                        /*

                          cout << "lambda1_2 " << lambda_1_lambda_2 << endl;
                          cout << "lambda1 " << lambda_1 << endl;
                               
                        */        
                        Vector<Rational> l1_nominator= centerpoint-vec0;
                        Vector<Rational> l1_denominator = vec2-vec0;
                        constraints.push_back({l1_nominator, l1_denominator, current, connection_index_current});
                        
                        /*
                          cout << "l1_nom: " << l1_nominator << endl;
                          cout << "l1_denom: " << l1_denominator << endl;
                        */

                        //Triangle 2: 
                       

                        //Computing Volume via determinant
                       
                        Rational lambda_2 = cross_product3(vec1-vec0, vec3-vec1)*(centerpoint-vec0);

                        Rational constant = cross_product3(vec1-vec0, vec1-vec0)*(centerpoint-vec0);

#if OPT_DEBUG
                        
                        
                          cout << "l2: " << lambda_2  << endl;
                          cout << "l1: " << constant << endl;

                          cout << qobj_coeff << endl;
                          cout << lobj_coeff << endl;
                          
                   
                        

                        

#endif

                         
                        if(lambda_2 + constant < 0)
                        {
                           
                            lobj_coeff[(lambda_index +1)%number_lambda] += Rational(1,6)*(-lambda_2);
                            cobj_coeff += Rational(1,6)*(-constant);
                            
                        }
                        else
                        {
                            lobj_coeff[(lambda_index +1)%number_lambda] += Rational(1,6)*(lambda_2);
                            cobj_coeff += Rational(1,6)*(constant);
                        }
                          

                        /*
                          cout << "lambda2 " << lambda_2 << endl;
                          cout << "const " << constant << endl;
                        */
                        
                        lambda_index++;


                        //___Determine the conter points to follow_____

                    
                        current_conter =  connection_index_current;
                        next_conter =  connection_index_next;
                    
                        first_conter = current_conter;
                        cout << "first conter (rectangle): " << current_conter << ", " << next_conter << endl;

                        //Get next conter element

                        long next_next_conter = next_conter;
                        
                        if(conter_boundary_array.size() == 2)
                        {
                            if(next_conter == conter_boundary_array[0])
                            {
                                next_next_conter = conter_boundary_array[1];
                            }
                            else
                            {
                                next_next_conter = conter_boundary_array[0];
                            }
                        }
                        else
                        { 
                            int64 next_conter_elements = 0;
                            for(int m = 0; m < conter_boundary_array.size(); m++)
                            {
                                if((total_adjacency(next_conter, conter_boundary_array[m])) && (conter_boundary_array[m] != current_conter))
                                {
                                    //NOTE: This path should bet taken only by one unique element!! No break statement!
                                    next_next_conter = conter_boundary_array[m];
                                    next_conter_elements++;
                                }    
                            }

                            //NOTE: NO next elment if subboundary has two elements...
                     
                            Assert(next_conter_elements <= 1);
                        }

                        current_conter = next_conter;
                        next_conter = next_next_conter;

                   
                    
                                    
                        
                    }


                }
               
                
                
            }
            else
            {
                cout << "current lambda after first step: " << lambda_index << endl;
                //Advance one step in subboundary
                long privious = current;
                current = next;
                int64 next_elements = 0;


                if(number_lambda == 5)
                {
                    int w = 0;
                }
                

               
                //Get next element
                if(dim == 2)
                {
                    next = privious;
                }
                else{
                    for(int l = 0; l < sub_boundary_array.size(); l++)
                    {
                        if((total_adjacency(current, sub_boundary_array[l])) && (sub_boundary_array[l] != privious))
                        {
                            //NOTE: This path should bet taken only by one unique element!! No break statement!
                            next = sub_boundary_array[l];
                            next_elements++;
                        }    
                    }

                    Assert(next_elements <= 1);

                }

                
                //NOTE: NO next elment if subboundary has two elements...
                //cout << "next: " << current << ", " << next << endl;
                cout << "sub: " << current << ", " << next << endl;
                
               

                
               
               
                bool32 walk_conter_boundary = true;
                bool32 walk_sub_boundary = false;
                bool32 cross_detected = false;
                bool32 skip_triangle = false;
                //NOTE: I intentionally set the conter values to sub values!
                int64 cross_current_conter = current;
                int64 cross_next_conter = next;

                    

                do
                {
                    //Walk on conter boundary and delay walk on sub boundary
                       
                        
                    //NOTE: The order of the cases matterss!

                    // Since in the 0ths step only the degenarte traingle case is checked (exept the subboundary is one point)
                    // We might have to do some work on the conterboundary first: see for example the standard 3 simplex
                    // The second condition on the if ensures that in the case where the conterboundary only contains one point, this case is skipped and the degenerate triangel case is taken
                    //The third and forth conditions make sure that
                    //1. a triangle in a cross is not walked through indefiniley
                    //2. the connecting lines (lambdas) get not repeated
                    
                    if( total_adjacency(current, next_conter) && (current_conter != next_conter)
                        && (current_conter != cross_current_conter) && (next_conter != cross_next_conter)
                        && (((current_conter != cross_next_conter) && (next_conter != cross_current_conter)) ||
                            (conter_points.size() != 2))
                        && !skip_triangle   )
                    {
                         /* case:
                   current |----------| current_conter
                           | \        |
                           |   \      |
                           |      \   |
                           |        \ |
                      next |----------| next_conter

                        */

                        //Cross!!
                        if(total_adjacency(next, current_conter))
                        {

                             //Register first cross
                            if(cross_detected == false)
                            {
                                cross_detected = true;
                                cross_current_conter = current_conter;
                                cross_next_conter = next_conter;
                            }


                            Vector<Rational> vecx1 = vertices.row(next) - vertices.row(current);
                            Vector<Rational> vecy2 = vertices.row(next_conter) -vertices.row(current);
                            Vector<Rational> vecz = centerpoint - vertices.row(current);

                            Rational cross_product_reverse_Z = cross_product3(vecy2, vecx1)*vecz;

                            Vector<Rational> vecy2Z = vertices.row(current_conter) -vertices.row(current);

                            Rational cross_product_Z = cross_product3(vecy2Z, vecx1)*vecz;


                            cout <<"reverse Z :" << convert_to<double>(cross_product_reverse_Z)<< endl;
                            cout << "Z :" << convert_to<double>(cross_product_Z) << endl;
                            if( cross_product_Z <= 0)
                            {
                                cout << convert_to<double>(cross_product_reverse_Z)<< endl;
                                cout << convert_to<double>(cross_product_Z) << endl;
                                int debug_W = 0;
                            }
                        
                            
                            //NOTE: If the first cross is the last surface we consider (in the subboundary)
                            //- and we run out of point in the conter boundary, so current_conter is first
                            // conter - The correct computation is a degenerat triangel not a triangle!
                            // afterwards we have to break
                            if(current_conter == first_conter)
                            {
                                //NOTE: if first surface is degenarate triangle we must be able to start,
                                // so no continue if k == 1
                                if(move_in_conter)
                                {
                                    //continue;
                                }

                                
                                if(move_in_conter){

                                if(cross_product_reverse_Z >= 0 && !reverse_dir)
                                {
                                    // continue;
                                }
                                if(cross_product_reverse_Z <= 0 && reverse_dir)
                                {
                                    //continue;
                                }
                                    }
                                
                                //NOTE: we remain at the same position in the conter_boundary
                                //but dont't get in this case anymore 
                             }


                            if((next_conter == first_conter) && !total_adjacency(next, next_conter))
                            {
                                //NOTE: THIS only works if the conterboundary is of size 2
                                //I let this like that, to see what happens in other cases
                                //continue;
                            }

                             if((next == first_sub) && !total_adjacency(next, next_conter))
                            {
                                
                                int64 next_next = getnextElement(current_conter, next_conter, conter_boundary_array, total_adjacency);

                                bool32 rectangle = false;

                                for(int l = 0; l < facets.rows(); l++)
                                {
                                    if(facets(l, current) && facets(l,next_next)
                                       && facets(l, next) && facets(l, next_conter) )
                                    {
                                        rectangle = true;
                                    }
                                }

                                
                                

                                
                                if(reverse_dir && !rectangle)
                                {
                                    ///skip_triangle = true;
                                
                                    //continue;
                                }

                               
                                
                            }

                             if(next == first_sub && next_conter == first_conter)
                             {
                                 //skip_triangle = true;
                                 //continue;
                             }

                            

                            

                           

                            

                            //NOTE: THE order vec2 and then vec1 is important.
                            //This must coinside with ter order in the computation of lambda1 in the (degen triangle) and rectangle case!
                             cout << convert_to<double>(cross_product_reverse_Z) << endl;
                             cout << convert_to<double>( vecx1) << endl;
                             cout << convert_to<double>(vecy2) << endl;
                             cout << convert_to<double>(vecz) << endl;

                             if(cross_product_Z > 0)
                            {

                                 bool32 is_facet = false;

                                for(int l = 0; l < facets.rows(); l++)
                                {
                                    if(facets(l, current) && facets(l,current_conter)
                                       && facets(l, next)  )
                                    {
                                        is_facet = true;
                                    }
                                }

                                if(is_facet)
                                {
                                    skip_triangle = true;
                                    continue;
                                }
                                
                            }
                            //else: we are in the correct path
                        }

                        //compute the volume and go to next_conter


                        Vector<Rational> vec0 = vertices.row(current);
                        Vector<Rational> vec1 = vertices.row(current_conter);
                        Vector<Rational> vec2 = vertices.row(next_conter);


                        //Computing Volume via determinant
                        Rational lambda_1_lambda_2 = cross_product3(vec1-vec0, vec2-vec0)*(centerpoint-vec0);

                        
                        qobj_coeff(lambda_index, (lambda_index +1)%number_lambda) = Rational(1,6)*AbsoluteValue(lambda_1_lambda_2);
                               
                                
                        Vector<Rational> l1_nominator= centerpoint-vec0;
                        Vector<Rational> l1_denominator = vec1-vec0;

                        constraints.push_back({l1_nominator, l1_denominator, current, current_conter});
                        lambda_index++;
                            
                            

                        cout << "conter (Triangle): " <<  current_conter << ", " << next_conter << endl;

                           

                        walk_conter_boundary = true;
                        walk_sub_boundary = false;
                        
                       
                    

                    }

                    else if(total_adjacency(current_conter, next))
                    {
                       
                         /* case:
                   current |----------| current_conter
                           |        / |
                           |      /   |
                           |   /      |
                           | /        |
                      next |----------| next_conter

                        */

                        //compute the volume and go to next
                        skip_triangle = false;

DEG_TRIAG:                  

                        //Compute Volume:
                    
                        //Triangle 1: 
                            Vector<Rational> vec0 = vertices.row(current);
                        Vector<Rational> vec1 = vertices.row(next);
                        Vector<Rational> vec2 = vertices.row(current_conter);


                        //Computing Volume via determinant
                        Rational lambda_1_lambda_2 = cross_product3(vec2-vec0, vec2-vec1)*(centerpoint-vec0);
                        Rational lambda_1 =cross_product3(vec2-vec0, vec1-vec0)*(centerpoint-vec0);
#if OPT_DEBUG
                        
                        
                          cout << "l1_l2: " << lambda_1_lambda_2  << endl;
                          cout << "l1: " << lambda_1 << endl;

                          cout << qobj_coeff << endl;
                          cout << lobj_coeff << endl;
                          
                   
                        

                        

#endif

                         
                        //NOTE in lobj_coeff += necessary since we might push linear terms ahead
                          
                        if(lambda_1_lambda_2 + lambda_1 < 0)
                        {
                            qobj_coeff(lambda_index, (lambda_index +1)%number_lambda) = Rational(1,6)*(-lambda_1_lambda_2);
                            lobj_coeff[lambda_index] += Rational(1,6)*(-lambda_1);
                            
                        }
                        else 
                        {
                            qobj_coeff(lambda_index, (lambda_index +1)%number_lambda) = Rational(1,6)*(lambda_1_lambda_2);
                            lobj_coeff[lambda_index] += Rational(1,6)*(lambda_1);
                        }
                        
                               
                                
                        Vector<Rational> l1_nominator= centerpoint-vec0;
                        Vector<Rational> l1_denominator = vec2-vec0;
                        constraints.push_back({l1_nominator, l1_denominator, current, current_conter});
                        
                        
                        //Triangle 2: 
                       

                        //Computing Volume via determinant
                       
                        Rational lambda_2 =cross_product3(vec1-vec0, vec2-vec1)*(centerpoint-vec0);
                        Rational constant = cross_product3(vec1-vec0, vec1-vec0)*(centerpoint-vec0);

#if OPT_DEBUG
                        
                        
                          cout << "l2: " << lambda_2  << endl;
                          cout << "l1: " << constant << endl;
                          
                          cout << qobj_coeff << endl;
                          cout << lobj_coeff << endl;
                        

                        

#endif
                          
      
                        if(lambda_2 + constant < 0)
                        {
                           
                            lobj_coeff[(lambda_index +1)%number_lambda] += Rational(1,6)*(-lambda_2);
                            cobj_coeff += Rational(1,6)*(-constant);
                            
                        }
                        else
                        {
                            lobj_coeff[(lambda_index +1)%number_lambda] += Rational(1,6)*(lambda_2);
                            cobj_coeff += Rational(1,6)*(constant);
                        }
                          

                       
                        
                        lambda_index++;
                    

                        //Delay the conter walk 
                   
                              cout << "conter (degenarate Triangle): " << current_conter << endl;
                        walk_conter_boundary = false;
                        walk_sub_boundary = true;
                            

                            
                    }

                    
                        
                    else 
                    {
                        /* case:
                   current |----------| current_conter
                           |          |
                           |          |
                           |          |
                           |          |
                      next |----------| next_conter

                        */

                        skip_triangle = false;

                        //compute volume and go to next and next_conter
                        Vector<Rational> vec0 = vertices.row(current);
                        Vector<Rational> vec1 = vertices.row(next);
                        Vector<Rational> vec2 = vertices.row(current_conter);
                        Vector<Rational> vec3 = vertices.row(next_conter);


                        //Triangle 1: 
                        //Computing Volume via determinant

                     
                        Rational lambda_1_lambda_2 = cross_product3(vec2-vec0, vec3-vec1)*(centerpoint-vec0);
                        Rational lambda_1 = cross_product3(vec2-vec0, vec1-vec0)*(centerpoint-vec0);

#if OPT_DEBUG
                        
                        
                          cout << "l1_l2: " << lambda_1_lambda_2  << endl;
                          cout << "l1: " << lambda_1 << endl;

                          cout << qobj_coeff << endl;
                          cout << lobj_coeff << endl;
                          
                   
                        

                        

#endif
                         
                        if(lambda_1_lambda_2 + lambda_1 < 0)
                        {
                            qobj_coeff(lambda_index, (lambda_index +1)%number_lambda) = Rational(1,6)*(-lambda_1_lambda_2);
                            lobj_coeff[lambda_index] += Rational(1,6)*(-lambda_1);
                            
                        }
                        else
                        {
                            qobj_coeff(lambda_index,(lambda_index +1)%number_lambda) = Rational(1,6)*(lambda_1_lambda_2);
                            lobj_coeff[lambda_index] += Rational(1,6)*(lambda_1);
                        }
                        

                                
                        Vector<Rational> l1_nominator= centerpoint-vec0;
                        Vector<Rational> l1_denominator = vec2-vec0;
                        constraints.push_back({l1_nominator, l1_denominator, current, current_conter});
                        

                        /*
                          cout << "vec0: " << vec0 << endl;
                          cout << "vec1: " << vec1 << endl;
                          cout << "vec2: " << vec2 << endl;
                          cout << "vec3: " << vec3 << endl;
                   
                          cout << "center: " << centerpoint << endl;

                          cout << "lambda1_2 " << lambda_1_lambda_2 << endl;
                          cout << "lambda1 " << lambda_1 << endl;
                               
                                
                    
                        

                          cout << "l1_nom: " << l1_nominator << endl;
                          cout << "l1_denom: " << l1_denominator << endl;

                        */

                        //Triangle 2: 
                       

                        //Computing Volume via determinant
                       
                        Rational lambda_2 = cross_product3(vec1-vec0, vec3-vec1)*(centerpoint-vec0);

                        Rational constant = cross_product3(vec1-vec0, vec1-vec0)*(centerpoint-vec0);

 #if OPT_DEBUG
                        
                        
                          cout << "l2: " << lambda_2  << endl;
                          cout << "l1: " << constant << endl;
                          cout << qobj_coeff << endl;
                          cout << lobj_coeff << endl;
                          
                   
                        

                        

#endif

                         
                        if(lambda_2 + constant < 0)
                        {
                           
                            lobj_coeff[(lambda_index +1)%number_lambda] += Rational(1,6)*(-lambda_2);
                            cobj_coeff += Rational(1,6)*(-constant);
                            
                        }
                        else
                        {
                            lobj_coeff[(lambda_index +1)%number_lambda] += Rational(1,6)*(lambda_2);
                            cobj_coeff += Rational(1,6)*(constant);
                        }
                          

                        /*
                          cout << "lambda2 " << lambda_2 << endl;
                          cout << "const " << constant << endl;
                        */
                       
                        
                        lambda_index++;

                            
                            
                        cout << "conter (rectangle): " << current_conter << ", " << next_conter << endl;

                        walk_conter_boundary = true;
                        walk_sub_boundary = true;
                    }



                   
                   

                    //Get next conter element
                    if(walk_conter_boundary == true)
                    {
                        long next_next_conter = next_conter;
                        
                        
                        if(conter_boundary_array.size() == 2)
                        {
                            if(next_conter == conter_boundary_array[0])
                            {
                                next_next_conter = conter_boundary_array[1];
                            }
                            else
                            {
                                next_next_conter = conter_boundary_array[0];
                            }
                        }
                        else
                        {
                            int64 next_conter_elements = 0;
                            for(int m = 0; m < conter_boundary_array.size(); m++)
                            {
                                if((total_adjacency(next_conter, conter_boundary_array[m])) && (conter_boundary_array[m] != current_conter))
                                {
                                    //NOTE: This path should bet taken only by one unique element!! No break statement!
                                    next_next_conter = conter_boundary_array[m];
                                    next_conter_elements++;
                                }    
                            }

                            //NOTE: NO next elment if subboundary has two elements...
                     
                            Assert(next_conter_elements <= 1);
                        }

                        move_in_conter = true;

                        current_conter = next_conter;
                        next_conter = next_next_conter;
                    }




               
                        

                        
                        
                }
                while(walk_conter_boundary && !walk_sub_boundary);



                if(next == first_sub)
                {

                    //Clean up the triangles at the first sub_boundary point
                   
                    while(current_conter != first_conter)
                    {
                
                        /* case:
                 first_sub |----------| current_conter
                           | \        |
                           |   \      |
                           |      \   |
                           |        \ |
                           |----------| next_conter

                        */

                        //compute the volume and go to next_conter


                        Vector<Rational> vec0 = vertices.row(first_sub);
                        Vector<Rational> vec1 = vertices.row(current_conter);
                        Vector<Rational> vec2 = vertices.row(next_conter);


                        //Computing Volume via determinant
                        Rational lambda_1_lambda_2 = cross_product3(vec1-vec0, vec2-vec0)*(centerpoint-vec0);
                        qobj_coeff(lambda_index, (lambda_index +1)%number_lambda) = Rational(1,6)*AbsoluteValue(lambda_1_lambda_2);
                               
                                
                        Vector<Rational> l1_nominator= centerpoint-vec0;
                        Vector<Rational> l1_denominator = vec1-vec0;

                        constraints.push_back({l1_nominator, l1_denominator, first_sub, current_conter});
                        lambda_index++;
                            
                            

                        cout << "conter (Triangle): " <<  current_conter << ", " << next_conter << endl;
                    

                        //Get next conter element
                        long next_next_conter = next_conter;
                        
                        
                        if(conter_boundary_array.size() == 2)
                        {
                            if(next_conter == conter_boundary_array[0])
                            {
                                next_next_conter = conter_boundary_array[1];
                            }
                            else
                            {
                                next_next_conter = conter_boundary_array[0];
                            }
                        }
                        else
                        {
                            int64 next_conter_elements = 0;
                            for(int m = 0; m < conter_boundary_array.size(); m++)
                            {
                                if((total_adjacency(next_conter, conter_boundary_array[m])) && (conter_boundary_array[m] != current_conter))
                                {
                                    //NOTE: This path should bet taken only by one unique element!! No break statement!
                                    next_next_conter = conter_boundary_array[m];
                                    next_conter_elements++;
                                }    
                            }

                            //NOTE: NO next elment if subboundary has two elements...
                     
                            Assert(next_conter_elements <= 1);
                        }

                        current_conter = next_conter;
                        next_conter = next_next_conter;

                        
                    }
                }        

      
               
            }
        }


        //Assert(constraints.size() == number_lambda);

















        Vector<double> print_centerpoint = convert_to<double>(centerpoint);
        cout << "(DEBUG) Centerpoint " << print_centerpoint << endl; 







        //NOTE: We can't use the rays as inequality here. We need to convert it to a h-representation!
       
        BigObject cone("Cone<Rational>");
        cone.take("INPUT_RAYS") << rays;
        Matrix<Rational> cone_h_representation = cone.give("FACETS");

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

        int scaling_exponent = static_cast<int>(highest_exponent) - 1;
        
        
        double tentimes = std::pow(10.0, scaling_exponent);
        
        if(tentimes > 0)
        {
            cone_h_representation /= convert_to<Rational>(tentimes);
        }




        //TODO: create model for QCQP solver

        
        bool32 clamp_precision = false;
        bool32 diffrent_contraint = false;
        Vector<Rational> first_direction = zero_vector<Rational>(3);
        Vector<Rational> second_direction = zero_vector<Rational>(3);
        double int_sol[number_lambda + dimension + 1];
        double    sol[number_lambda + dimension];
        double    objval;
        double int_objval;
            
        Assert(constraints.size() == number_lambda);

#if OPT_DEBUG            
        cout << "The model objective is: " << endl;
        cout << convert_to<double>(qobj_coeff) << endl;
        cout << convert_to<double>(lobj_coeff) << endl;
        cout << convert_to<double>(cobj_coeff) << endl;
#endif


        //gurobi test:

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
           
            
        //NOTE: The bounds on the lambda are very important otherwise the computation takes forever... should i add the norm constraint for u to?

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
        //lambda
        error = GRBaddvars(model, number_lambda, 0, NULL, NULL, NULL, obj_lin, lambda_lb, lambda_ub, NULL , NULL);
        
        
        
        if (error) goto QUIT1;

        //u
        error = GRBaddvars(model, dimension, 0, NULL, NULL, NULL, NULL, u_lb, u_ub, NULL, NULL);
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
            
            
            
            
 

        /* Hyperplane constraint for each lambda */

        for(int32 w = 0; w < constraints.size(); w++)
        {
            for(int32 v = 0; v < dimension; v++)
            {
                //Note number_lambda + v is the u indicex
                qrow[v] = w; qcol[v] = number_lambda+v; qval[v] = convert_to<double>(constraints[w].denum_vec[v]);

                lind[v] = number_lambda+v; lval[v] = - convert_to<double>(constraints[w].num_vec[v]);
            }

#if OPT_DEBUG
                
            cout << "lambd for " << constraints[w].start_vec_index << " " << constraints[w].end_vec_index <<endl;
            cout << "num =" <<  convert_to<double>(constraints[w].num_vec) << endl;
            cout << "denum =" <<  convert_to<double>(constraints[w].denum_vec) << endl;
#endif
                
            error = GRBaddqconstr(model, dimension, lind, lval, dimension, qrow, qcol, qval, GRB_EQUAL, 0, NULL);
            // error = GRBaddqconstr(model, dimension, lind, lval, dimension, qrow, qcol, qval, GRB_LESS_EQUAL, 0.1, NULL);
            //error = GRBaddqconstr(model, dimension, lind, lval, dimension, qrow, qcol, qval, GRB_GREATER_EQUAL, -0.1, NULL);
            if (error) goto QUIT1;

                                     
        }
            
               
#if OPT_DEBUG            
        cout << convert_to<double>(rays) << endl;
        cout << convert_to<double>(cone_h_representation) << endl;
#endif
            

        for(int32 w = 0; w < cone_h_representation.rows(); w++)
        {

            for(int32 v = 0; v < dimension; v++)
            {
                //Note number_lambda + v is the u index
                if(clamp_precision == true)
                {
                    lind[v] = number_lambda+v; lval[v] = decimal_truncate(convert_to<double>(cone_h_representation(w,v)),3);
                    /*
                    cout << "lind " << v << " " << lval[v] << " " << convert_to<double>(cone_h_representation(w,v)) << endl;
                    real64 truncate_test = truncate_to_digits_v2(convert_to<double>(cone_h_representation(w,v)),3);
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

           

        
        if(diffrent_contraint == false)
        {
            //Norm constraint for u, or just u > 0
            //u2_qrow[number_lambda] = 1.0;
            

            for(int32 w = 0; w < dimension; w++)
            {
                //Note number_lambda + w is the u index
                u2_qrow[w] = number_lambda+w; u2_qcol[w] = number_lambda+w; u2_qval[w]= 1.0;
            }
            
            error = GRBaddqconstr(model, 0, NULL, NULL, dimension, u2_qrow, u2_qcol, u2_qval, GRB_EQUAL, 1.0, "norm");
            if (error) goto QUIT1;
        }
        else
        {
            /*

            char ztype[dimension];
            for (int32 w = 0; w < dimension; w++)
            {
                ztype[w] = GRB_BINARY;
            }
            error = GRBaddvars(model, dimension, 0, NULL, NULL, NULL, NULL,  NULL, NULL, ztype, NULL);
            if (error) goto QUIT1;

            error = GRBupdatemodel(model);
            if (error) goto QUIT1;

            // Add big-M constraints: -M*z[i] <= x[i] <= M*z[i]
            for (int32 w = 0; w < dimension; w++)
            {
                int ind[2];
                double val[2];
                double M = 1.0;

                // x[i] - M*z[i] <= 0
                ind[0] = number_lambda + w;         // x[i]
                ind[1] = number_lambda + dimension + w;    // z[i]
                val[0] = 1.0;
                val[1] = -M;
                error = GRBaddconstr(model, 2, ind, val, GRB_LESS_EQUAL, 0.0, NULL);
                if (error) goto QUIT1;

                // x[i] + M*z[i] >= 0  <=> -x[i] - M*z[i] <= 0
                val[0] = -1.0;
                val[1] = -M;
                error = GRBaddconstr(model, 2, ind, val, GRB_LESS_EQUAL, 0.0, NULL);
                if (error) goto QUIT1;
            }

            int ind[dimension];
            double val[dimension];
            for (int32 w = 0; w < dimension; w++) {
                ind[w] = number_lambda + dimension + w;  // z[i]
                val[w] = 1.0;
            }

            error = GRBaddconstr(model, dimension, ind, val, GRB_GREATER_EQUAL, 1.0, "nonzero");
            if (error) goto QUIT1;
            */

           
            error = GRBaddvar(model, 0, NULL, NULL, 0.0, 0.5, GRB_INFINITY, GRB_CONTINUOUS, "aux_var");
            if (error) goto QUIT1;

            error = GRBupdatemodel(model);
            if (error) goto QUIT1;


            /* aux_var = 1-norm(u1, u2, u3) */
            int ind[] = {number_lambda, number_lambda + 1, number_lambda + 2};
            error = GRBaddgenconstrNorm(model, "1-norm", number_lambda + dimension, 3, ind, 1.0);
            if (error) goto QUIT1;

        }

            
       
            

        /* Optimize model */

        error = GRBoptimize(model);
        if (error) goto QUIT1;

           

        /* Capture solution information */

        

        if(diffrent_contraint == false)
        {
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

                first_direction = minimal_directions.row(i);

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

        }

        else
        {
            error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
            if (error) goto QUIT1;

            error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &int_objval);
            if (error) goto QUIT1;

            error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, number_lambda + dimension + 1, int_sol);
            if (error) goto QUIT1;

            printf("\nOptimization complete\n");
            if (optimstatus == GRB_OPTIMAL) {
                printf("Optimal objective: %.4e\n", objval);

                for(int32 w = 0; w < number_lambda; w++)
                {
                    cout << "l" << w << "=" << int_sol[w]  << " ";
                   
                }

                //NOTE: WE cast the direction back to Rational; Do I have to worry about precision??
                
                for(int32 w = 0; w < dimension; w++)
                {
                    second_direction[w] = convert_to<Rational>(int_sol[number_lambda + w]);
                    cout << "u" << w << "=" << int_sol[number_lambda + w]  << " ";
                }

                cout << std::defaultfloat << std::setprecision(6);
                cout << endl;

                cout << "First Direction: " << convert_to<double>(first_direction) << endl;
                cout << "Second Direction: " << convert_to<double>(second_direction) << endl;
                
                
                
            } else if (optimstatus == GRB_INF_OR_UNBD) {
                printf("Model is infeasible or unbounded\n");
            } else {
                printf("Optimization was stopped early\n");
            }


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
            }
            
            exit(1);
        }

        

        /* Free model */

        GRBfreemodel(model);

        /* Free environment */

        GRBfreeenv(env);

        // For quadratic constraint

        Matrix<double> real_vertices = convert_to<double>(vertices);
        Matrix<double> real_conter_vertices = convert_to<double>(conter_vertices);
            
        for(int l = 0; l < constraints.size(); l++)
        {
                
            Vector<double> intersection_vec = real_vertices.row(constraints[l].start_vec_index) + sol[l]*(real_vertices.row(constraints[l].end_vec_index) - real_vertices.row(constraints[l].start_vec_index));
            real_conter_vertices /= intersection_vec;

#if OPT_DEBUG            
            cout << intersection_vec << endl;
#endif
                
                
        }

        //NOTE: It seems that volume computation is only possible for Rational vertices

        Matrix<Rational> rational_conter_vertices = convert_to<Rational>(real_conter_vertices);

        BigObject conter_polytope("Polytope<Rational>");

        conter_polytope.take("VERTICES") << (ones_vector<Rational>() | rational_conter_vertices);
             
        real64 ConterVolume = conter_polytope.give("VOLUME");

#if OPT_DEBUG           
        cout << "Real conter vertices: " << real_conter_vertices << endl;
        cout << "sub_vertices: " << sub_vertices << endl;
#endif

            
        real64 sub_volume = objval + FixVolume;
           
        minimal_volume_signature_indices.push_back(i);
        minimal_volumes.push_back(sub_volume);
           
            
        real64 total_volume = ConterVolume + objval + FixVolume;
            
        cout << "Volume: " << ConterVolume << "(conter) + " << objval << "(objval) + " << FixVolume << "(fix-sub) = " << total_volume << endl;


        // For norm constraint

        Matrix<double> int_real_vertices = convert_to<double>(vertices);
        Matrix<double> int_real_conter_vertices = convert_to<double>(conter_vertices);
            
        for(int l = 0; l < constraints.size(); l++)
        {
                
            Vector<double> int_intersection_vec = int_real_vertices.row(constraints[l].start_vec_index) + int_sol[l]*(int_real_vertices.row(constraints[l].end_vec_index) -int_real_vertices.row(constraints[l].start_vec_index));
            int_real_conter_vertices /= int_intersection_vec;

#if OPT_DEBUG            
            cout << int_intersection_vec << endl;
#endif
                
                
        }

        //NOTE: It seems that volume computation is only possible for Rational vertices

        Matrix<Rational> int_rational_conter_vertices = convert_to<Rational>(int_real_conter_vertices);

        BigObject int_conter_polytope("Polytope<Rational>");

        int_conter_polytope.take("VERTICES") << (ones_vector<Rational>() | int_rational_conter_vertices);
             
        real64 int_ConterVolume = int_conter_polytope.give("VOLUME");

#if OPT_DEBUG           
        cout << "Real conter vertices: " << int_real_conter_vertices << endl;
        cout << "sub_vertices: " << sub_vertices << endl;
#endif

            
       
            
        real64 int_total_volume = int_ConterVolume + int_objval + FixVolume;
            
        cout << "int_Volume: " << int_ConterVolume << "(conter) + " << int_objval << "(objval) + " << FixVolume << "(fix-sub) = " << int_total_volume << endl;


        

//NOTE: THIS is DEBUG code

        /*
        if((total_volume > 1001) || (total_volume < 999))
        {
            cout << "THE COMPUTED VOLUMES ARE NOT CORRECT !!!!";
            Assert(1 == 0);
                   
        }

        */
/*       
        if((total_volume > 167) || (total_volume < 166))
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



Vector<Rational> interior_points_methode(BigObject* p, Matrix<Rational> vertices, Vector<Rational> point,
                                         Matrix<Rational> Inequalities, memory_arena *Arena)
{
    
    halfspace start_halfspace = minimum_halfspace(p, vertices, point);
    Vector<Rational> direction = - start_halfspace.direction;
    real64 optimal_value = start_halfspace.value;
    real64 volume_difference = start_halfspace.value;

    interior_point *first_point = PushStruct(Arena, interior_point);
    first_point->point_x = convert_to<real64>(point[0]);
    first_point->point_y = convert_to<real64>(point[1]);
    first_point->point_z = convert_to<real64>(point[2]);
    first_point->direction_x = convert_to<real64>(direction[0]);
    first_point->direction_y = convert_to<real64>(direction[1]);
    first_point->direction_z = convert_to<real64>(direction[2]);
    first_point->value = optimal_value;
    first_point->volume_difference = volume_difference;
   
    
    cout << "Direction " << convert_to<double>(direction) << endl;

    cout << "First point: " << point << " ------ Halfspace-depth: " << optimal_value << endl;

    //Check for feasible
    //NOTE: The direction vector is noramlized!
    int64 outer_loop_iteration = 0;

    while( (outer_loop_iteration < 8))
    {
        int64 accuracy = 1;
        Rational lambda =  Rational(1, accuracy);
        Vector<Rational> start_point = point;
        

        while(!is_feasible(&Inequalities, start_point + lambda*direction))
        {
            accuracy++;
            lambda = Rational(1, accuracy);
        }

        Vector<Rational> first_test_point = start_point + lambda*direction;

        
        halfspace test_halfspace = minimum_halfspace(p, vertices, first_test_point);
        if(test_halfspace.value > optimal_value)
        {
            point = first_test_point;
            volume_difference = test_halfspace.value - optimal_value;
            direction = - test_halfspace.direction;
            optimal_value = test_halfspace.value;

            

            interior_point *log_point = PushStruct(Arena, interior_point);
            log_point->point_x = convert_to<real64>(point[0]);
            log_point->point_y = convert_to<real64>(point[1]);
            log_point->point_z = convert_to<real64>(point[2]);
            log_point->direction_x = convert_to<real64>(direction[0]);
            log_point->direction_y = convert_to<real64>(direction[1]);
            log_point->direction_z = convert_to<real64>(direction[2]);
            log_point->value = optimal_value;
            log_point->volume_difference = volume_difference;
            log_point->step_size = 1;
            log_point->path = 1;
            
           
        }
            
        else
        {
        
            int64 inner_loop_iterations = 0;
            int64 number_steps = 100;

           

            
            Rational step_size = lambda/Rational(number_steps,1);
            Vector<Rational> step_search_point = start_point;
            Vector<Rational> step_search_direction = direction;
            real64 step_search_optimal_value = optimal_value;
            int64 found_index;

            bool32 found = false;

            
            interior_point *log_point = PushStruct(Arena, interior_point);
            log_point->step_size = 100;


            

           
            for(int j = 0; j < number_steps; j++)
            {            

                //NOTE: All this points shoud be feasible
                Vector<Rational> test_point = start_point + (j*step_size)*direction;
                Assert(is_feasible(&Inequalities, test_point));

                //halfspace step_test_halfspace = minimum_halfspace(p, vertices, test_point);

                /*
                cout << step_test_halfspace.value << endl;
                if(step_test_halfspace.value > 70.0)
                {
                    cout <<  convert_to<double>(test_point) << endl;
                    cout << step_size << endl;
                    cout << outer_loop_iteration << endl;
                    cout <<  convert_to<double>(direction) << endl;

                    halfspace ext_halfspace = minimum_halfspace(p, vertices, test_point);

                    cout << ext_halfspace.value << endl;
                    
                    Vector<double> test = convert_to<double>(test_point);

                    halfspace bad_halfspace = minimum_halfspace(p, vertices, convert_to<Rational>(test));

                    cout << bad_halfspace.value << endl;
                    cout << convert_to<double>(convert_to<Rational>(test)) << endl;

                    cout << test_point << endl;
                    cout << convert_to<Rational>(test) << endl;


                }
                */
                Vector<double> d_test_point = convert_to<double>(test_point);
                
                Vector<Rational> rounded_test_point = zero_vector<Rational>(3);

                /*
                rounded_test_point[0] = convert_to<Rational>(round_to_digits(d_test_point[0],9));
                rounded_test_point[1] = convert_to<Rational>(round_to_digits(d_test_point[1],9));
                rounded_test_point[2] = convert_to<Rational>(round_to_digits(d_test_point[2],9));
                */
                rounded_test_point[0] = convert_to<Rational>(d_test_point[0]);
                rounded_test_point[1] =  convert_to<Rational>(d_test_point[1]);
                rounded_test_point[2] =  convert_to<Rational>(d_test_point[2]);

                halfspace d_step_test_halfspace =  minimum_halfspace(p, vertices, rounded_test_point);

                log_point->line_values[j] = d_step_test_halfspace.value;
                log_point->line_points[j].X = convert_to<double>(rounded_test_point[0]);
                log_point->line_points[j].Y = convert_to<double>(rounded_test_point[1]);
                log_point->line_points[j].Z = convert_to<double>(rounded_test_point[2]);
                log_point->line_directions[j].X = convert_to<real64>(d_step_test_halfspace.direction[0]);
                log_point->line_directions[j].Y = convert_to<real64>(d_step_test_halfspace.direction[1]);
                log_point->line_directions[j].Z = convert_to<real64>(d_step_test_halfspace.direction[2]);
               
                
                if(d_step_test_halfspace.value > step_search_optimal_value + 0.001)
                {

                    //Avoid badly condition search_directions
                    //OR randomly choose better one???
                  
                    step_search_optimal_value = d_step_test_halfspace.value;
                    step_search_direction = - d_step_test_halfspace.direction;
                    step_search_point = rounded_test_point;
                    found = true;
                    found_index = j;
                    log_point->found = true;
                    log_point->path = 2;
                  
                   
                
                }
                else
                {
                    //NOTE: Optimization: the objective function is convex
                    
                    if(found)
                    {
                        //break;
                    }
                     if((j == 1) && (step_search_optimal_value + 0.001 > d_step_test_halfspace.value))
                        {
                            step_search_optimal_value = d_step_test_halfspace.value;
                            step_search_direction = - convert_to<double>(d_step_test_halfspace.direction);
                            step_search_point = test_point;
                            found_index = j;
                            found = true;
                            
                            log_point->path = 6;
                        }

                    
                }
                
                inner_loop_iterations++;
            
            }

            
           
                    
            point = step_search_point;
            volume_difference = step_search_optimal_value - optimal_value;
            direction = step_search_direction;
            optimal_value = step_search_optimal_value;

            //smooth the subgradient to avoid unreasonalble small step_sizes
            Vector<Rational> smooth_direction = zero_vector<Rational>(3);
            real64  smooth_direction_x = 0;
            real64  smooth_direction_y = 0;
            real64  smooth_direction_z = 0;

            int64 start_index = 0;
            if(found_index - 7 >= 0)
            {
                start_index = found_index - 7;
            }

            for(int j= 0; j < 15; j++)
            {
                
                smooth_direction_x += (-log_point->line_directions[start_index + j].X);
                smooth_direction_y += (-log_point->line_directions[start_index + j].Y);
                smooth_direction_z += (-log_point->line_directions[start_index + j].Z);      
            }

            smooth_direction_x = smooth_direction_x/(15.0);
            smooth_direction_y = smooth_direction_y/(15.0);
            smooth_direction_z = smooth_direction_z/(15.0);

            smooth_direction[0] = convert_to<Rational>(smooth_direction_x);
            smooth_direction[1] = convert_to<Rational>(smooth_direction_y);
            smooth_direction[2] = convert_to<Rational>(smooth_direction_z);

            direction = smooth_direction;

            

            log_point->point_x = convert_to<real64>(point[0]);
            log_point->point_y = convert_to<real64>(point[1]);
            log_point->point_z = convert_to<real64>(point[2]);
            log_point->direction_x = convert_to<real64>(direction[0]);
            log_point->direction_y = convert_to<real64>(direction[1]);
            log_point->direction_z = convert_to<real64>(direction[2]);
            log_point->value = optimal_value;
            log_point->volume_difference = volume_difference;
          
           
  
            

        }
        
        cout << "Iteration " << outer_loop_iteration << ": " << "Point: " << point  << " ------ Halfspace-depth: " << optimal_value << endl;  

        outer_loop_iteration++;
    }

    Vector<double> centerpoint = convert_to<double>(point);

    cout << "RESULT (Iterations " << outer_loop_iteration << "): -- (approx.) Centerpoint: " << centerpoint << "    -- halfspace depth: " << optimal_value << endl;

    cout << "---------- SUMMARY -----------" << endl;


    std::vector<real64> exact_checker;
     std::vector<real64> checker;

     for(int i = 0; i < (Arena->Used/sizeof(interior_point)); i++)
    {
        interior_point* step = (interior_point *) (Arena->Base + i*sizeof(interior_point));
       
        Vector<double> test_exact = zero_vector<double>(3);
        Vector<double> test = zero_vector<double>(3);

        real64 xx = step->point_x;
        real64 yy = step->point_y;
        real64 zz = step->point_z;

        test_exact[0] = xx;
        test_exact[1] = yy;
        test_exact[2] = zz;

        cout << test_exact << endl;

        real64 x = round_to_digits(step->point_x, 9);
        real64 y = round_to_digits(step->point_y, 9);
        real64 z = round_to_digits(step->point_z, 9);


        test[0] = x;
        test[1] = y;
        test[2] = z;

        cout << test << endl;


        halfspace exact_test_halfspace =  minimum_halfspace(p, vertices, convert_to<Rational>(test_exact));

        exact_checker.push_back(exact_test_halfspace.value);

        halfspace test_halfspace =  minimum_halfspace(p, vertices, convert_to<Rational>(test));

        checker.push_back(test_halfspace.value);
          
    }

    for(int i = 0; i < (Arena->Used/sizeof(interior_point)); i++)
    {
        interior_point* step = (interior_point *) (Arena->Base + i*sizeof(interior_point));
        cout << "Iteration " << i << ": Point: " << step->point_x << " " << step->point_y << " " << step->point_z << "  -----> (next) direction: " << step->direction_x << " " << step->direction_y << " " << step->direction_z << " (value: " << step->value << " dif: " << step->volume_difference << " step_size: " << step->step_size << "with: " << step->path << ") "<<endl;

        

        cout <<"DEBUG checksum = " << exact_checker.at(i) << " ( " << checker.at(i) << " )" << endl;

        if(step->step_size == 100)
        {
            for(int j = 0; j < 100; j++)
            {
                cout << "            (" << step->line_points[j].X << " " << step->line_points[j].Y << " " <<  step->line_points[j].Z << ")  -----" << step->line_values[j] << "     (" << step->line_directions[j].X << " " << step->line_directions[j].Y << " " <<  step->line_directions[j].Z << ")" << endl;  
            }

        }
       
       
          
    }
    
    return point;    
}






Vector<Rational> interior_points_methode_v2(BigObject* p, Matrix<Rational> vertices, Vector<Rational> point,
                                         Matrix<Rational> Inequalities, memory_arena *Arena)
{
    
    halfspace start_halfspace = minimum_halfspace(p, vertices, point);
    Vector<double> direction = - convert_to<double>(start_halfspace.direction);
    real64 optimal_value = start_halfspace.value;
    real64 volume_difference = start_halfspace.value;

    Vector<double> d_point = convert_to<double>(point);

    interior_point *first_point = PushStruct(Arena, interior_point);
    first_point->point_x = d_point[0];
    first_point->point_y = d_point[1];
    first_point->point_z = d_point[2];
    first_point->direction_x = direction[0];
    first_point->direction_y =direction[1];
    first_point->direction_z = direction[2];
    first_point->value = optimal_value;
    first_point->volume_difference = volume_difference;
   
    
    cout << "Direction " << direction << endl;

    cout << "First point: " << point << " ------ Halfspace-depth: " << optimal_value << endl;

    //Check for feasible
    //NOTE: The direction vector is noramlized!
    int64 outer_loop_iteration = 0;

    while( (outer_loop_iteration < 8))
    {
        int64 accuracy = 1;
        double lambda =  1;
        Vector<double> start_point = d_point;
        

        while(!is_feasible(&Inequalities, convert_to<Rational>(start_point + lambda*direction)))
        {
            accuracy++;
            lambda = 1 / (real64) accuracy;
        }

        Vector<double> first_test_point = start_point + lambda*direction;

        
        halfspace test_halfspace = minimum_halfspace(p, vertices, convert_to<Rational>(first_test_point));
        if(test_halfspace.value > optimal_value)
        {
            d_point = first_test_point;
            volume_difference = test_halfspace.value - optimal_value;
            direction = - convert_to<double>(test_halfspace.direction);
            optimal_value = test_halfspace.value;

            

            interior_point *log_point = PushStruct(Arena, interior_point);
            log_point->point_x = d_point[0];
            log_point->point_y = d_point[1];
            log_point->point_z = d_point[2];
            log_point->direction_x = direction[0];
            log_point->direction_y = direction[1];
            log_point->direction_z = direction[2];
            log_point->value = optimal_value;
            log_point->volume_difference = volume_difference;
            log_point->step_size = 1;
            log_point->path = 1;
            
           
        }
            
        else
        {
        
            int64 inner_loop_iterations = 0;
            int64 number_steps = 100;

           

            
            double step_size = lambda/ (double)number_steps;
            Vector<double> step_search_point = start_point;
            Vector<double> step_search_direction = direction;
            real64 step_search_optimal_value = optimal_value;
            int64 found_index;

            bool32 found = false;

            
            interior_point *log_point = PushStruct(Arena, interior_point);
            log_point->step_size = 100;


            

           
            for(int j = 0; j < number_steps; j++)
            {            

                //NOTE: All this points shoud be feasible
                Vector<double> test_point = start_point + (j*step_size)*direction;
                Assert(is_feasible(&Inequalities, convert_to<Rational>(test_point)));

               

                halfspace d_step_test_halfspace =  minimum_halfspace(p, vertices, convert_to<Rational>(test_point));

                log_point->line_values[j] = d_step_test_halfspace.value;
                log_point->line_points[j].X = test_point[0];
                log_point->line_points[j].Y = test_point[1];
                log_point->line_points[j].Z = test_point[2];
                log_point->line_directions[j].X = convert_to<real64>(d_step_test_halfspace.direction[0]);
                log_point->line_directions[j].Y = convert_to<real64>(d_step_test_halfspace.direction[1]);
                log_point->line_directions[j].Z = convert_to<real64>(d_step_test_halfspace.direction[2]);
               
                
                if(d_step_test_halfspace.value > step_search_optimal_value + 0.001)
                {

                    //Avoid badly condition search_directions
                    //OR randomly choose better one???
                  
                    step_search_optimal_value = d_step_test_halfspace.value;
                    step_search_direction = - convert_to<double>(d_step_test_halfspace.direction);
                    step_search_point = test_point;
                    found = true;
                    found_index = j;
                    log_point->found = true;
                    log_point->path = 2;
                  
                   
                
                }
                else
                {
                    //NOTE: Optimization: the objective function is convex
                    
                    if(found)
                    {
                        //break;
                    }
                    
                     if((j == 1) && (step_search_optimal_value + 0.001 > d_step_test_halfspace.value))
                        {
                            step_search_optimal_value = d_step_test_halfspace.value;
                            step_search_direction = - convert_to<double>(d_step_test_halfspace.direction);
                            step_search_point = test_point;
                            found_index = j;
                            found = true;
                            
                            log_point->path = 6;
                        }

                    
                }
                
                inner_loop_iterations++;
            
            }

             
           
                    
            d_point = step_search_point;
            volume_difference = step_search_optimal_value - optimal_value;
            direction = step_search_direction;
            optimal_value = step_search_optimal_value;

            //smooth the subgradient to avoid unreasonalble small step_sizes
            Vector<double> smooth_direction = zero_vector<double>(3);
            real64  smooth_direction_x = 0;
            real64  smooth_direction_y = 0;
            real64  smooth_direction_z = 0;

            int64 start_index = 0;
            if(found_index - 7 >= 0)
            {
                start_index = found_index - 7;
            }

            for(int j= 0; j < 15; j++)
            {
                
                smooth_direction_x += (-log_point->line_directions[start_index + j].X);
                smooth_direction_y += (-log_point->line_directions[start_index + j].Y);
                smooth_direction_z += (-log_point->line_directions[start_index + j].Z);      
            }

            smooth_direction_x = smooth_direction_x/(15.0);
            smooth_direction_y = smooth_direction_y/(15.0);
            smooth_direction_z = smooth_direction_z/(15.0);

            smooth_direction[0] = smooth_direction_x;
            smooth_direction[1] = smooth_direction_y;
            smooth_direction[2] = smooth_direction_z;

            direction = smooth_direction;

            

            log_point->point_x = d_point[0];
            log_point->point_y = d_point[1];
            log_point->point_z = d_point[2];
            log_point->direction_x = direction[0];
            log_point->direction_y = direction[1];
            log_point->direction_z = direction[2];
            log_point->value = optimal_value;
            log_point->volume_difference = volume_difference;
          
           
  
            

        }
        
        cout << "Iteration " << outer_loop_iteration << ": " << "Point: " << point  << " ------ Halfspace-depth: " << optimal_value << endl;  

        outer_loop_iteration++;
    }

    Vector<double> centerpoint = convert_to<double>(point);

    cout << "RESULT (Iterations " << outer_loop_iteration << "): -- (approx.) Centerpoint: " << centerpoint << "    -- halfspace depth: " << optimal_value << endl;

    cout << "---------- SUMMARY -----------" << endl;


    std::vector<real64> exact_checker;
     std::vector<real64> checker;

     for(int i = 0; i < (Arena->Used/sizeof(interior_point)); i++)
    {
        interior_point* step = (interior_point *) (Arena->Base + i*sizeof(interior_point));
       
        Vector<double> test_exact = zero_vector<double>(3);
        Vector<double> test = zero_vector<double>(3);

        real64 xx = step->point_x;
        real64 yy = step->point_y;
        real64 zz = step->point_z;

        test_exact[0] = xx;
        test_exact[1] = yy;
        test_exact[2] = zz;

        cout << test_exact << endl;

        real64 x = round_to_digits(step->point_x, 9);
        real64 y = round_to_digits(step->point_y, 9);
        real64 z = round_to_digits(step->point_z, 9);


        test[0] = x;
        test[1] = y;
        test[2] = z;

        cout << test << endl;


        halfspace exact_test_halfspace =  minimum_halfspace(p, vertices, convert_to<Rational>(test_exact));

        exact_checker.push_back(exact_test_halfspace.value);

        halfspace test_halfspace =  minimum_halfspace(p, vertices, convert_to<Rational>(test));

        checker.push_back(test_halfspace.value);
          
    }

    for(int i = 0; i < (Arena->Used/sizeof(interior_point)); i++)
    {
        interior_point* step = (interior_point *) (Arena->Base + i*sizeof(interior_point));
        cout << "Iteration " << i << ": Point: " << step->point_x << " " << step->point_y << " " << step->point_z << "  -----> (next) direction: " << step->direction_x << " " << step->direction_y << " " << step->direction_z << " (value: " << step->value << " dif: " << step->volume_difference << " step_size: " << step->step_size << "with: " << step->path << ") "<<endl;

        

        cout <<"DEBUG checksum = " << exact_checker.at(i) << " ( " << checker.at(i) << " )" << endl;

        if(step->step_size == 100)
        {
            for(int j = 0; j < 100; j++)
            {
                cout << "            (" << step->line_points[j].X << " " << step->line_points[j].Y << " " <<  step->line_points[j].Z << ")  -----" << step->line_values[j] << "     (" << step->line_directions[j].X << " " << step->line_directions[j].Y << " " <<  step->line_directions[j].Z << ")" << endl;  
            }

        }
       
       
          
    }
    
    return point;    
}



Vector<Rational>sub_grad_descend(BigObject* p, Matrix<Rational> vertices, Vector<Rational> point,
                                         Matrix<Rational> Inequalities, memory_arena *Arena)
{

    int64 Num_Iteration = 20;
    Vector<Rational> current_point = point;
    halfspace start_halfspace = minimum_halfspace(p, vertices, current_point);
    Vector<Rational> current_direction = - start_halfspace.direction;
    int64 unfeasable_count = 0;
    
   

    log_point *first_point = PushStruct(Arena, log_point);
    first_point->point.X = convert_to<real64>(current_point[0]);
    first_point->point.Y = convert_to<real64>(current_point[1]);
    first_point->point.Z = convert_to<real64>(current_point[2]);
    first_point->direction.X = convert_to<real64>(current_direction[0]);
    first_point->direction.Y = convert_to<real64>(current_direction[1]);
    first_point->direction.Z = convert_to<real64>(current_direction[2]);
    first_point->value = start_halfspace.value;

    
    for (int i = 1; i < Num_Iteration + 1; i++)
    {
        Vector<Rational> test_point = current_point + Rational(1,i)*current_direction;

        //Check for feasible
        //NOTE: The direction vector is noramlized!
        int64 outer_loop_iteration = 0;

        if(!is_feasible(&Inequalities, test_point))
        {
            unfeasable_count++;

            int64 accuracy = 1;
            Rational lambda =  Rational(1, accuracy);
      
            while(!is_feasible(&Inequalities, test_point))
            {
                accuracy++;
                lambda = Rational(1, accuracy);
                test_point = current_point + lambda*Rational(1,i)*current_direction;
            }


            Assert(is_feasible(&Inequalities, test_point));

        }

        /*
        //NOTE: This procedure was implemented before the scaling of the cone-h-representation
        // The new code is in the beginning of the minimum_halfspace_slices function
        // It helped in some cases... but has no logical reason to work
          Vector<double> d_test_point = convert_to<double>(test_point);
                
          Vector<Rational> rounded_test_point = zero_vector<Rational>(3);

        
          rounded_test_point[0] = convert_to<Rational>(round_to_digits(d_test_point[0],9));
          rounded_test_point[1] = convert_to<Rational>(round_to_digits(d_test_point[1],9));
          rounded_test_point[2] = convert_to<Rational>(round_to_digits(d_test_point[2],9));
        
          rounded_test_point[0] = convert_to<Rational>(d_test_point[0]);
          rounded_test_point[1] =  convert_to<Rational>(d_test_point[1]);
          rounded_test_point[2] =  convert_to<Rational>(d_test_point[2]);

          current_point = rounded_test_point;
        */

        current_point = test_point;
        
        halfspace test_halfspace = minimum_halfspace(p, vertices, current_point);
        current_direction = - test_halfspace.direction;
        
       
    
   

        log_point *intermediate_point = PushStruct(Arena, log_point);
        intermediate_point->point.X = convert_to<real64>(current_point[0]);
        intermediate_point->point.Y = convert_to<real64>(current_point[1]);
        intermediate_point->point.Z = convert_to<real64>(current_point[2]);
        intermediate_point->direction.X = convert_to<real64>(current_direction[0]);
        intermediate_point->direction.Y = convert_to<real64>(current_direction[1]);
        intermediate_point->direction.Z = convert_to<real64>(current_direction[2]);
        intermediate_point->value = test_halfspace.value;

    }


    for(int i = 0; i < (Arena->Used/sizeof(log_point)); i++)
    {
        log_point* step = (log_point *) (Arena->Base + i*sizeof(log_point));
        cout << "Iteration " << i << ": Point: " << step->point.X << " " << step->point.Y
             << " " << step->point.Z << "  -----> (next) direction: " << step->direction.X << " " << step->direction.Y << " " << step->direction.Z << " (value: " << step->value << ") "<<endl;
     
    }

    cout << "(infeasabilities " << unfeasable_count << ")" << endl;

    return current_point;

}


//TOOLs


Matrix<Rational> scale(double xfactor, double yfactor, double zfactor)
{
    Matrix<Rational> scaling_matrix = zero_matrix<Rational>(3,3);
    scaling_matrix(0,0) = convert_to<Rational>(xfactor);
    scaling_matrix(1,1) = convert_to<Rational>(yfactor);
    scaling_matrix(2,2) = convert_to<Rational>(zfactor);

    return scaling_matrix;
    
}


Matrix<Rational> rotate_y(double angle)
{
        
    Matrix<Rational> rotation_matrix_x = zero_matrix<Rational>(3,3);
    rotation_matrix_x(0,0) = convert_to<Rational>(cos(angle));
    rotation_matrix_x(0,2) = convert_to<Rational>(sin(angle));
    rotation_matrix_x(1,1) = 1;
    rotation_matrix_x(2,0) = -convert_to<Rational>(sin(angle));
    rotation_matrix_x(2,2) = convert_to<Rational>(cos(angle));

    return rotation_matrix_x;
}

Matrix<Rational> rotate_z(double angle)
{
        
    Matrix<Rational> rotation_matrix_x = zero_matrix<Rational>(3,3);
    rotation_matrix_x(0,0) = convert_to<Rational>(cos(angle));
    rotation_matrix_x(0,1) = convert_to<Rational>(sin(angle));
    rotation_matrix_x(1,0) = - convert_to<Rational>(sin(angle));
    rotation_matrix_x(1,1) = convert_to<Rational>(cos(angle));
    rotation_matrix_x(2,2) = 1;

    return rotation_matrix_x;
}


inline double norm(Vector<Rational> vec)
{
    return sqrt(convert_to<double>(sqr(vec)));
}
 


int main(int argc, const char* argv[])
{
    void* BaseAddress = (void *)(0);
    uint64  TotalSize = Megabytes(1);

    memory_arena Arena = getMemoryArena(BaseAddress, TotalSize);
    
    Main pm;
    pm.set_application("fan");
    
   

    //3D
    Matrix<Rational> simplex3_vertices_old = {{0,  0,  0},{10, 0,  0},{0, 10, 0},{0, 0, 10}};

    Matrix<Rational> cube3zero_vertices_old = {{5,  5,  5},{5, -5,  5},{5,  5, -5},{5, -5, -5},
                                           {-5,  5,  5},{-5, -5,  5},{-5,  5, -5},{-5, -5, -5}};
    Matrix<Rational> cube3_vertices_old = {{10,  10,  10},{10, 0, 10},{10,  10, 0},{10, 0, 0},
                                           {0,  10,  10},{0, 0, 10},{0,  10, 0},{0, 0, 0}};
    Matrix<Rational> twopyramid3_vertices_old = {{5, 5, 0},{5, -5,  0},{-5, 5, 0},{-5, -5, 0},
                                                 {0, 0, 5},{0,0,-5}};


    Array<Matrix<Rational>> slices_simplex3_vertices(3);

    slices_simplex3_vertices[0] =  {{0, 0, 0}, {0, 10, 0}, {0, 0, 10}};
    slices_simplex3_vertices[1] =  {{1, 0, 0}, {1, 9, 0}, {1, 0, 9}};
    slices_simplex3_vertices[2] =  {{2, 0, 0}, {2, 8, 0}, {2, 0, 8}};
    
    /*
    slices_simplex3_vertices[3] =  {{3, 0, 0}, {3, 7, 0}, {3, 0, 7}};
    slices_simplex3_vertices[4] =  {{4, 0, 0}, {4, 6, 0}, {4, 0, 6}};
    slices_simplex3_vertices[5] =  {{5, 0, 0}, {5, 5, 0}, {5, 0, 5}};
    slices_simplex3_vertices[6] =  {{6, 0, 0}, {6, 4, 0}, {6, 0, 4}};
    slices_simplex3_vertices[7] =  {{7, 0, 0}, {7, 3, 0}, {7, 0, 3}};
    slices_simplex3_vertices[8] =  {{8, 0, 0}, {8, 2, 0}, {8, 0, 2}};
    slices_simplex3_vertices[9] =  {{9, 0, 0}, {9, 1, 0}, {9, 0, 1}};
    slices_simplex3_vertices[10] =  {{10, 0, 0}};
    */

    Array<Matrix<Rational>> slices_cube3_vertices(3);

    slices_cube3_vertices[0] =  {{0, 0, 0}, {0, 10, 0}, {0, 10, 10}, {0, 0, 10}};
    slices_cube3_vertices[1] =  {{1, 0, 0}, {1, 10, 0}, {1, 10, 10}, {1, 0, 10}};
    slices_cube3_vertices[2] =  {{2, 0, 0}, {2, 10, 0}, {2, 10, 10}, {2, 0, 10}};

    /*
    slices_cube3_vertices[3] =  {{3, 0, 0}, {3, 10, 0}, {3, 10, 10}, {3, 0, 10}};
    slices_cube3_vertices[4] =  {{4, 0, 0}, {4, 10, 0}, {4, 10, 10}, {4, 0, 10}};
    slices_cube3_vertices[5] =  {{5, 0, 0}, {5, 10, 0}, {5, 10, 10}, {5, 0, 10}};
    slices_cube3_vertices[6] =  {{6, 0, 0}, {6, 10, 0}, {6, 10, 10}, {6, 0, 10}};
    slices_cube3_vertices[7] =  {{7, 0, 0}, {7, 10, 0}, {7, 10, 10}, {7, 0, 10}};
    slices_cube3_vertices[8] =  {{8, 0, 0}, {8, 10, 0}, {8, 10, 10}, {8, 0, 10}};
    slices_cube3_vertices[9] =  {{9, 0, 0}, {9, 10, 0}, {9, 10, 10}, {9, 0, 10}};
    slices_cube3_vertices[10] = {{10, 0, 0}, {10, 10, 0}, {10, 10, 10}, {10, 0, 10}};
    */
  
       
    Array<Matrix<Rational>> simplex3_vertices(1);
    simplex3_vertices[0] = simplex3_vertices_old;

    Array<Matrix<Rational>> cube3zero_vertices(1);
    cube3zero_vertices[0] = cube3zero_vertices_old;

    Array<Matrix<Rational>> cube3_vertices(1);
    cube3_vertices[0] = cube3_vertices_old;

    Array<Matrix<Rational>> twopyramid3_vertices(1);
    twopyramid3_vertices[0] = twopyramid3_vertices_old;



    //3D Degenerate Imput
    //Only one 2d slice in 3d --> should be given in 2d!

    Array<Matrix<Rational>> one_2slices_simplex3_vertices(1);

    one_2slices_simplex3_vertices[0] =  {{0, 0, 0}, {0, 10, 0}, {0, 0, 10}};

    Array<Matrix<Rational>> one_2slices_cube3_vertices(1);

    one_2slices_cube3_vertices[0] =  {{0, 0, 0}, {0, 10, 0}, {0, 10, 10}, {0, 0, 10}};
    


    //4D
    //Array<Matrix<Rational>> simplex4(1);

    //simplex4[0] = {{0, 0, 0, 0}, {10, 0, 0,0}, {0,10,0,0},{0, 0, 10,0}, {0, 0,0, 10}};
    
    //Vector<Rational> cube_centerpoint = {1, 1, 1,1};
    //sub_grad_descend_slices(simplex4, cube_centerpoint, 0,  &Arena);

   

    //Custom Test
    Array<Matrix<Rational>> oert_slices_simplex3_vertices(3);

    oert_slices_simplex3_vertices[0] =  {{0, 0, 0}, {0, 10, 0}, {0, 0, 10}};
    oert_slices_simplex3_vertices[1] =  {{1, 0, 0}, {1, 9, 0}, {1, 0, 9},{1,9,9}};
    oert_slices_simplex3_vertices[2] =  {{2, 0, 0}, {2, 8, 0}, {2, 0, 8}};

    
    //3d:
    //TEST
    std::mt19937 gen(76);                         // Fixed seed
    std::uniform_real_distribution<double> dist1(0.0, 10.0);
    std::uniform_real_distribution<double> dist2(0.0, 10.0);
    std::uniform_real_distribution<double> dist3(0.0, 10.0);
     std::uniform_real_distribution<double> dist4(0.0, 10.0);  
     // std::uniform_real_distribution<double> dist5(0.0, 1.0);
     //std::uniform_real_distribution<double> dist6(0.0, 7.0);
    // std::uniform_real_distribution<double> dist7(0.0, 1.0);
     
    

    Array<Matrix<Rational>> slices(4);
    //slice 1

    for(int i = 0; i < 4; i++)
    {
        if(i == 0)
        {
            Rational slice_1_1_y = convert_to<Rational>(dist1(gen));
            Rational slice_1_1_z = convert_to<Rational>(dist1(gen));

            Vector<Rational> slice_1_1= zero_vector<Rational>(3);
            slice_1_1[0] = i;
            slice_1_1[1] = slice_1_1_y;
            slice_1_1[2] = slice_1_1_z;

            Rational slice_1_2_y = convert_to<Rational>(dist1(gen));
            Rational slice_1_2_z = convert_to<Rational>(dist1(gen));

            Vector<Rational> slice_1_2= zero_vector<Rational>(3);
            slice_1_2[0] = i;
            slice_1_2[1] = slice_1_2_y;
            slice_1_2[2] = slice_1_2_z;

       
            Rational slice_1_3_y = convert_to<Rational>(dist1(gen));
            Rational slice_1_3_z = convert_to<Rational>(dist1(gen));
        
            Vector<Rational> slice_1_3= zero_vector<Rational>(3);
            slice_1_3[0] = i;
            slice_1_3[1] = slice_1_3_y;
            slice_1_3[2] = slice_1_3_z;

            Matrix<Rational> slice_1 = vector2row(slice_1_1) / slice_1_2 / slice_1_3;
        
            slices[i] = slice_1;
        }
        if(i == 1)
        {
            Rational slice_1_1_y = convert_to<Rational>(dist2(gen));
            Rational slice_1_1_z = convert_to<Rational>(dist2(gen));

            Vector<Rational> slice_1_1= zero_vector<Rational>(3);
            slice_1_1[0] = i;
            slice_1_1[1] = slice_1_1_y;
            slice_1_1[2] = slice_1_1_z;

            Rational slice_1_2_y = convert_to<Rational>(dist2(gen));
            Rational slice_1_2_z = convert_to<Rational>(dist2(gen));

            Vector<Rational> slice_1_2= zero_vector<Rational>(3);
            slice_1_2[0] = i;
            slice_1_2[1] = slice_1_2_y;
            slice_1_2[2] = slice_1_2_z;

       
            Rational slice_1_3_y = convert_to<Rational>(dist2(gen));
            Rational slice_1_3_z = convert_to<Rational>(dist2(gen));
        
            Vector<Rational> slice_1_3= zero_vector<Rational>(3);
            slice_1_3[0] = i;
            slice_1_3[1] = slice_1_3_y;
            slice_1_3[2] = slice_1_3_z;

            Matrix<Rational> slice_1 = vector2row(slice_1_1) / slice_1_2 / slice_1_3;
        
            slices[i] = slice_1;

        }
         if(i == 2)
        {
            Rational slice_1_1_y = convert_to<Rational>(dist3(gen));
            Rational slice_1_1_z = convert_to<Rational>(dist3(gen));

            Vector<Rational> slice_1_1= zero_vector<Rational>(3);
            slice_1_1[0] = i;
            slice_1_1[1] = slice_1_1_y;
            slice_1_1[2] = slice_1_1_z;

            Rational slice_1_2_y = convert_to<Rational>(dist3(gen));
            Rational slice_1_2_z = convert_to<Rational>(dist3(gen));

            Vector<Rational> slice_1_2= zero_vector<Rational>(3);
            slice_1_2[0] = i;
            slice_1_2[1] = slice_1_2_y;
            slice_1_2[2] = slice_1_2_z;

       
            Rational slice_1_3_y = convert_to<Rational>(dist3(gen));
            Rational slice_1_3_z = convert_to<Rational>(dist3(gen));
        
            Vector<Rational> slice_1_3= zero_vector<Rational>(3);
            slice_1_3[0] = i;
            slice_1_3[1] = slice_1_3_y;
            slice_1_3[2] = slice_1_3_z;

            Matrix<Rational> slice_1 = vector2row(slice_1_1) / slice_1_2 / slice_1_3;
        
            slices[i] = slice_1;

        }
         
          if(i == 3)
        {
            Rational slice_1_1_y = convert_to<Rational>(dist4(gen));
            Rational slice_1_1_z = convert_to<Rational>(dist4(gen));

            Vector<Rational> slice_1_1= zero_vector<Rational>(3);
            slice_1_1[0] = i;
            slice_1_1[1] = slice_1_1_y;
            slice_1_1[2] = slice_1_1_z;

            Rational slice_1_2_y = convert_to<Rational>(dist4(gen));
            Rational slice_1_2_z = convert_to<Rational>(dist4(gen));

            Vector<Rational> slice_1_2= zero_vector<Rational>(3);
            slice_1_2[0] = i;
            slice_1_2[1] = slice_1_2_y;
            slice_1_2[2] = slice_1_2_z;

       
            Rational slice_1_3_y = convert_to<Rational>(dist4(gen));
            Rational slice_1_3_z = convert_to<Rational>(dist4(gen));
        
            Vector<Rational> slice_1_3= zero_vector<Rational>(3);
            slice_1_3[0] = i;
            slice_1_3[1] = slice_1_3_y;
            slice_1_3[2] = slice_1_3_z;

            Matrix<Rational> slice_1 = vector2row(slice_1_1) / slice_1_2 / slice_1_3;
        
            slices[i] = slice_1;

        }
          /*
          if(i == 4)
        {
            Rational slice_1_1_y = convert_to<Rational>(dist5(gen));
            Rational slice_1_1_z = convert_to<Rational>(dist5(gen));

            Vector<Rational> slice_1_1= zero_vector<Rational>(3);
            slice_1_1[0] = i;
            slice_1_1[1] = slice_1_1_y;
            slice_1_1[2] = slice_1_1_z;

            Rational slice_1_2_y = convert_to<Rational>(dist5(gen));
            Rational slice_1_2_z = convert_to<Rational>(dist5(gen));

            Vector<Rational> slice_1_2= zero_vector<Rational>(3);
            slice_1_2[0] = i;
            slice_1_2[1] = slice_1_2_y;
            slice_1_2[2] = slice_1_2_z;

       
            Rational slice_1_3_y = convert_to<Rational>(dist5(gen));
            Rational slice_1_3_z = convert_to<Rational>(dist5(gen));
        
            Vector<Rational> slice_1_3= zero_vector<Rational>(3);
            slice_1_3[0] = i;
            slice_1_3[1] = slice_1_3_y;
            slice_1_3[2] = slice_1_3_z;

            Matrix<Rational> slice_1 = vector2row(slice_1_1) / slice_1_2 / slice_1_3;
        
            slices[i] = slice_1;

        }
         */
             
    
    

    }

    


    Array<BigObject> hyperplanes(10);
    
    for(int i = 0; i < 10; i++)
    {
        Matrix<Rational> h_ineq = {{i,-1,0,0},{-i,1,0,0}};
  
        BigObject h("Polytope<Rational>");
        h.take("INEQUALITIES") << h_ineq;
        hyperplanes[i] = h;

    }

    //Matrix<Rational> test2 =  vector2row(slices[2].row(0)) / slices[2].row(2);

    /*
    slices[2].row(0)[1] =  convert_to<Rational>(2.29577);
    slices[2].row(0)[2] =  convert_to<Rational>(4.43453);

    slices[2].row(1)[1] =  convert_to<Rational>(9.13962);
    slices[2].row(1)[2] =  convert_to<Rational>(5.34414);

    slices[2].row(2)[1] =  convert_to<Rational>(4.30699);
    slices[2].row(2)[2] =  convert_to<Rational>(4.57205);

    
    */

    Matrix<Rational> all_slices = slices[0];
   
    
    for(int i = 1; i < 4; i++)
    {
        all_slices /= slices[i];
    }
    /*
    
     Matrix<Rational> test_slices =  {{0, 0, 0}, {0, 10, 0}, {0, 0, 10},
                   {1, 0, 0}, {1, 9, 0}, {1, 0, 9},
                   {2, 0, 0}, {2, 8, 0}, {2, 0, 8},
                   {3, 0, 0}, {3, 7, 0}, {3, 0, 7},
                   {4, 0, 0}, {4, 6, 0}, {4, 0, 6},
                   {5, 0, 0}, {5, 5, 0}, {5, 0, 5},
                   {6, 0, 0}, {6, 4, 0}, {6, 0, 4},
                   {7, 0, 0}, {7, 3, 0}, {7, 0, 3},
                   {8, 0, 0}, {8, 2, 0}, {8, 0, 2},
                   {9, 0, 0}, {9, 1, 0}, {9, 0, 1}
                   };*/


     // 90: 1.5708
     // 60: 1.0472
     // 45: 0.7854
     // 30: 0.5236
    
    
     // Matrix<Rational> rotation_z = rotate_z(  1.5708);
     // Matrix<Rational> rotation = rotate_z(0.7508);
    /*
     Matrix<Rational> full_rotation =T(rotation_z);

    Vector<Rational> cc = {1,3,3};
    
     all_slices = test_slices * full_rotation;
     Vector<Rational> cc_test = cc *T(rotation_z);
   
   

    cout << convert_to<double>(all_slices) << endl;
    cout << "rotateted Centerpoint " << convert_to<double>(cc_test) << endl;

    */

    //Array<Matrix<Rational>> random_slices_vertices(3);
    //Array<BigObject> convex_body_h(3);


    BigObject convex_body("Polytope<Rational>");
    convex_body.take("POINTS") <<  (ones_vector<Rational>() | all_slices);

    int num_slices = 0;
    for(int i = 0; i < 4; i++)
    {
        BigObject convex_body_h1 = call_function("intersection", convex_body, hyperplanes[i]);
        if(convex_body_h1.give("FEASIBLE"))
        {
            num_slices++;     
        }
    }

    Array<Matrix<Rational>> random_slices_vertices(num_slices);
    Array<BigObject> convex_body_h(num_slices);


    
    for(int i = 0; i < 4; i++)
    {
        BigObject convex_body_h1 = call_function("intersection", convex_body, hyperplanes[i]);
        if(convex_body_h1.give("FEASIBLE"))
        {
        Matrix<Rational> matrix_h1 = convex_body_h1.give("VERTICES");
        random_slices_vertices[i] = matrix_h1.minor(All, range_from(1));
        convex_body_h[i] = convex_body_h1;
        }
    }
  
    int ss = 1;
    
    Vector<Rational> point = convex_body_h[ss].give("REL_INT_POINT");
    Vector<Rational> new_centerpoint = point.slice(range_from(1));

    
    cout << "slice 1 : " << endl;
    cout << convert_to<double>(random_slices_vertices[0]) << endl;

    cout << "slice 2: " << endl;
    cout << convert_to<double>(random_slices_vertices[1]) << endl;

    //cout << "slice 3 : " << endl;
    //cout << convert_to<double>(random_slices_vertices[2]) << endl;

    

    cout << "centerpoint : " << convert_to<double>(new_centerpoint) << endl;


    /*
    new_centerpoint[1] = convert_to<Rational>(2.10);
    new_centerpoint[2] = convert_to<Rational>(4.8);
    */

    
    Array<Matrix<Rational>> full_test(1);
    Matrix<Rational> full_matrix = convex_body.give("VERTICES");
    full_test[0] = full_matrix.minor(All,range_from(1));
    Vector<Rational> full_point =  convex_body.give("REL_INT_POINT");;
    Vector<Rational> new_full_centerpoint = full_point.slice(range_from(1));
    
    
    cout << "Convex body points: " << endl;

    //cout << convert_to<double>(full_test[0]) << endl;

    //convex_body.call_method("VISUAL");




    /*
    
    log_point full_best = sub_grad_descend_slices(full_test, new_full_centerpoint, 0,  &Arena);

     Vector<Rational> full_best_centerpoint = zero_vector<Rational>(3);
    full_best_centerpoint[0] = convert_to<Rational>(full_best.point.X);
    full_best_centerpoint[1] = convert_to<Rational>(full_best.point.Y);
    full_best_centerpoint[2] = convert_to<Rational>(full_best.point.Z);

    cout << "Best_Centerpoint" << convert_to<double>(full_best_centerpoint) << endl;
    cout << endl;

    Assert(1 == 0);
    */
    
    

 
    


    
    
    /*
    new_full_centerpoint[0] = convert_to<Rational>(1);
    new_full_centerpoint[1] = convert_to<Rational>(4);
    new_full_centerpoint[2] = convert_to<Rational>(5);
    */
    
   
    log_point best = sub_grad_descend_slices(random_slices_vertices, new_centerpoint, ss,  &Arena);
    

    Vector<Rational> best_centerpoint = zero_vector<Rational>(3);
    best_centerpoint[0] = convert_to<Rational>(best.point.X);
    best_centerpoint[1] = convert_to<Rational>(best.point.Y);
    best_centerpoint[2] = convert_to<Rational>(best.point.Z);

    cout << "Best_Centerpoint" << convert_to<double>(best_centerpoint) << endl;
    cout << endl;
    
    //Combinatroics:
    cout << "---COMBINATORICS---" << endl;
    int n_vertices = convex_body.give("N_VERTICES");
    cout << "Number of Vertices: " << n_vertices << endl;
    Vector<Integer> f_vector = convex_body.give("F_VECTOR");
    cout << "f-vector: " << f_vector << endl;
    bool32 simple = convex_body.give("SIMPLE");
    cout << "simple: " << simple << endl;
    bool32 simplicial = convex_body.give("SIMPLICIAL");
    cout << "simplical: " << simplicial << endl;


    cout << endl;
    cout << "---GEOMETRY---" << endl;
    //Geometry:
    
    Vector<Rational> hom_centroid = convex_body.give("CENTROID");
    Vector<Rational> centroid = hom_centroid.slice(range_from(1));

    
    cout << "Centroid:    " << convert_to<double>(centroid) << endl;
    cout << "Norm( Centroid - C ):    " << norm(centroid-best_centerpoint) << endl;
    cout << "Proj_Norm (Centroid -C):    " <<  norm(centroid.slice(range_from(1))-best_centerpoint.slice(range_from(1))) << endl;
 
    Vector<Rational> hom_steiner_point = convex_body.give("STEINER_POINT");
    Vector<Rational> steiner_point = hom_steiner_point.slice(range_from(1));
    cout << "Steiner Point " << convert_to<double>(steiner_point) << endl;
    cout << "Norm( Steiner - C ):    " << norm(steiner_point-best_centerpoint) << endl;
    

    std::pair<Rational,Vector<Rational>> minimal_ball_pair = convex_body.give("MINIMAL_BALL");
    Vector<Rational> hom_minimal_ball_center = minimal_ball_pair.second;
    Vector<Rational> minimal_ball_center = hom_minimal_ball_center.slice(range_from(1));
    cout << "Minimal-Ball Center:    " << convert_to<double>(minimal_ball_center) << endl;
    cout << "Norm( MB-Center - C ):    " << norm(minimal_ball_center-best_centerpoint) << endl;
    

    Vector<Rational> hom_vertex_barycenter = convex_body.give("VERTEX_BARYCENTER");
    Vector<Rational> vertex_barycenter = hom_vertex_barycenter.slice(range_from(1));
    cout << "Vertex_barycenter:    " << convert_to<double>(vertex_barycenter) << endl;
    cout << "Norm( V-Bary - C ):   " << norm(vertex_barycenter-best_centerpoint) << endl;
    real32 min_vert_angle =  convex_body.give("MINIMAL_VERTEX_ANGLE");
    cout << "minimal vertex_angle:    " << min_vert_angle << endl;

    
      
    //TEST regular polytope
    //Don't have to rotate them!!!
    // cross, cube, 24 cell(4d) are the only ones with rational points
     
    //BigObject unit_cube = call_function("cube", 3);
    //BigObject cube = call_function("scale", unit_cube, 10);

    // cout << cube.give("VERTICES") << endl;
    //rotate  0, 30, 45; -> 6 times -> translate

    // BigObject cross = call_function("cross", 5);

    //Can i do it with the rational quadratic extension??
    //BigObject dodecahedron = call_function("dodecahedron");

    

    //Test general poytopes 
          //simplex ...
          //cyclic(Int d, Int n)

          //pyramid with base as n-gon
          //rand_box

          //zonotope
   
    //rand_box     

     

  
    
    
   

    

  
 


    //Assert(1 == 0);

    

    
    
    //Matrix<Rational> rotation_matrix_x = {{1, 0, 0},}
    

    //BigObject p("Polytope<Rational>");
   
    //TODO: DO POINT INSTAED if vertices slices
    // p.take("VERTICES") << (ones_vector<Rational>() | cube3_vertices_old);

   
   

    //NOTE: for the interior point test point 3,3,3
    //
    //Vector<Rational> centerpoint = {0, 0, 2};
    //Vector<Rational> centerpoint = {0, 1, 2};

    /*
    centerpoint[0] = Rational(12,5);
    centerpoint[1] = Rational(-9,3);
    centerpoint[2] = Rational(-9,4);
    */
    
    /*
    centerpoint[0] = convert_to<Rational>(2.39269);
    centerpoint[1] = convert_to<Rational>(0.57302);
    centerpoint[2] = convert_to<Rational>(1.05122);
    */
      
    //Matrix<Rational> Inequalities = p.give("FACETS");

    //real64 Volume = p.give("VOLUME");
    //cout << "Volume: " << Volume << endl;

    /*

    if(!is_feasible(&Inequalities, centerpoint))
    {
        cout << "Point " << centerpoint << " is not contained in convex body" << endl;
        Assert(false);
    }
    

    cout << "Centerpoint" <<  centerpoint << endl;
    */

    //centerpoint = {2,2,3};
    
    //sub_grad_descend_slices(cube3_vertices, centerpoint, 0,  &Arena);

    

    // minimum_halfspace(&p, vertices, centerpoint);
    //interior_points_methode(&p, vertices, centerpoint, Inequalities, &Arena);



    //sub_grad_descend(&p, cube3_vertices_old, centerpoint, Inequalities, &Arena);
    
    //sub_grad_descend_slices(cube3_vertices, centerpoint, 0,  &Arena);

    //Testing degenerate input.

    //sub_grad_descend_slices(one_2slices_cube3_vertices, centerpoint, 0,  &Arena);

    //Vector<Rational> centerpoint_v2 = {1, 2, 2};
    //minimum_halfspace_v2(&p, vertices, centerpoint);
    //minimum_halfspace_slices(vertices, centerpoint, 11);
    
    //!!!!!
    //centerpoint = {1, 4, 1};
    //sub_grad_descend_slices(oert_slices_simplex3_vertices, centerpoint, 1,  &Arena);
    

    /*
    double test = 1;
    double test2 = 2;;
    int sol_size = 3;

    read_baron_results(&test, &test2, sol_size);
    */
    // minimum_halfspace(&p, vertices, centerpoint_v2);

}


