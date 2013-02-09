#include <Eigen/Dense>
#include <time.h>
#include "ndarray/eigen.h"



namespace gptest {

//using namespace Eigen;

template <typename dty>
double EuclideanDistance(dty*,dty*,int);

template <typename dtyi, typename dtyo>
dtyo ExpCovariogram(dtyi*,dtyi*,int);

template <typename datatype>
void merge_sort(datatype*,int*,int);

template <typename datatype>
int merge_scanner(datatype*,int*,int,int);

template <typename datatype>
class kd{
  private:
    int **tree,pts,dim,room,roomstep,*inn,masterparent;
    int *neighcand,neighfound,neighwant;
    datatype **data,*tosort;
    double *ddcand;
    
    void organize(int*,int,int,int);
    int find_node(datatype*);
    void nn_explore(datatype*,int,int);
    void walk_up(int,int,int);
    
   double (*distance)(datatype*,datatype*,int);
    
  public:
    ~kd();
    kd(int,int,datatype**,double(*)(datatype*,datatype*,int));
    void nn_srch(datatype*,int,int*,double*);

    void add_pt(datatype*);
    
    void get_tree(int,int*);
    void test_sort();
    void test_scanner();
    void black_box_test();
    

  
};


template <typename dtyi, typename dtyo>
class gaussianprocess{

  private:
    int pts,n_nn,maxmin,dim,room,called_interp,*neigh,calleddummy;

    double *ddneigh;
    dtyo *fn,*ggq,kriging_parameter,*lambda;
    dtyi **data,*max,*min,*vv;
    
    double **ggl,**gglin;
    
    Eigen::Matrix <dtyo,Eigen::Dynamic,Eigen::Dynamic> gg,ggin;
    kd<dtyi> *kptr;
          
    double (*distance)(dtyi*,dtyi*,int);
    dtyo (*covariogram)(dtyi*,dtyi*,int);
      
  public:
  
     double etime,ltime,ect,lct;
    
  
     ~gaussianprocess();
      
     gaussianprocess(float);
     
     gaussianprocess(int,int,ndarray::Array<dtyi,2,2>,ndarray::Array<dtyo,1,1>);
      
     /*gaussianprocess(int,int,\
     ndarray::Array<dtyi,2,2>,ndarray::Array<dtyi,2,2>,ndarray::Array<dtyi,2,2>,dtyo*);
     */
     
     //note: code will remember whether or not you input with maxs and
     //mins
     
     dtyo interpolate(ndarray::Array<dtyi,1,1>,ndarray::Array<dtyo,1,1>,int);
    
     void set_kp(int);
     void print_nn(ndarray::Array<int,1,1>);
     void set_lambda(dtyo);
     void print_ggrow(int,ndarray::Array<dtyo,1,1>);

};

}
