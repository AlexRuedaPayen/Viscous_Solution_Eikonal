#include "EF2d-base.hpp"
#include <fstream>
#include <utility>

Mesh2d::Mesh2d(const char * filename)
{
  std::ifstream  f(filename); 
  assert( f); 
  int unused, I[4] ; 
  f >> nv >> nt >> unused ;
  assert( f.good());
  t = new Simplex[nt];
  v = new Vertex[nv];
  assert( t && v); 
  double mes =0; 
  for(int i=0;i<nv;++i)
    { 
      f >> v[i] ; 
      assert( f.good());
    }

  for(int k=0;k<nt;++k)
    { 
      for(int i=0;i< 4; ++i)
	f >> I[i] ; 
      assert( f.good());
      t[k].build(v,I,-1);
      mes += t[k].mes; 
    }
  std::cout<< " End read " << nv << " " << nt << " mes =" << mes << std::endl; 
}

std::pair<int,int> wrap_edge(int a) const {
	int m1,m2;
	if (a%3==0) {
		m1=this->operator()(a/3,1);
		m2=this->operator()(a/3,2);
	}
	if (a%3==1) {
		m1=this->operator()(a/3,0);
		m2=this->operator()(a/3,2);
	}
	if (a%3==2) {
		m1=this->operator()(a/3,0);
		m2=this->operator()(a/3,1);
	}
	return(std::make_pair(std::min(m1,m2),std::max(m1,m2)));
}


void Initialization() {

	std::pair<int,int> fbuffer;
	working_edges=new Edges[2*nt+1];

	for (int i=0;i<3*nt;++i) {
	
		fbuffer=wrap_edge(i); 

		Edges * it=registered_edges[fbuffer.first];
		Edges * current_working_edge=(t[i/3].adjacent_edges[i%3]);							

		if (it==nullptr) {
			registered_edges[fbuffer.first]=current_working_edge;
		}

		while(it!=nullptr && (it->vertices).second<fbuffer.second) {
			it=it->next_get_edge;
		}

		//if another neighbour element has already been browsed


		if (it==nullptr) {
			registered_edges[fbuffer.first]=current_working_edge;
		}

		if ((it->vertices).second==fbuffer.second) {

			Edges * current_working_edge=(t[i/3].adjacent_edges[i%3]);

			current_working_edge->next_get_shell=it;
			image[findex].to_itself_fill_next=it;
		}


		if (it!=nullptr && (it->vertices).second==fbuffer.second) {

		
			register_edges[findex].to_preimage=preimage+x;	
			image[findex].to_itself_extern_next=begin_extern;
			image[findex].to_itself_fill_next=it;
			
			else {
				if (it!=nullptr) {
					it=it->to_itself_fill_prev;
					image[findex].to_itself_fill_prev=it;
			}
				if (image[findex].to_itself_fill_prev!=nullptr) {
				(image[findex].to_itself_fill_prev)->to_itself_fill_next=image+findex;
			}
			}	
			if (image[findex].to_itself_fill_next!=nullptr) {
			(image[findex].to_itself_fill_next)->to_itself_fill_prev=image+findex;
		}
			findex+=1;
		}
}
