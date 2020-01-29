#include "EF2d-base.hpp"
#include <fstream>
#include <utility>
#include <limits>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <set>

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

  initialized=0; 
}

std::pair<int,int> Mesh2d::wrap_edge(int a) const {
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


void Mesh2d::initialization() {
	register_edges=new Edges*[nv];
	for (int i=0;i<nv;++i) {
		register_edges[i]=nullptr;
	}
	working_edges=new Edges[3*nt]; //majoration extrême dans un cadre connexe on a seulement au max 2*nt+1

	int findex=0;
	std::pair<int,int> fbuffer;
	Edges * it1;
	Edges * it2;

	for (int i=0;i<3*nt;++i) {

		fbuffer=wrap_edge(i);
		it1=register_edges[fbuffer.first];
		if (it1==nullptr) {
			register_edges[fbuffer.first]=working_edges+findex;
			working_edges[findex].vertices=fbuffer;
  			working_edges[findex].first_neighbour=(i/3);
			findex+=1;
			continue;
		}

		if ((it1->vertices).second>fbuffer.second) {
			register_edges[fbuffer.first]=working_edges+findex;
			working_edges[findex].vertices=fbuffer;
			working_edges[findex].first_neighbour=(i/3);

			working_edges[findex].next_edge=it1;
			findex+=1;
			continue;
		}

		while(it1!=nullptr && (it1->vertices).second<fbuffer.second) {
			it2=it1;
			it1=it1->next_edge;
		}

		if (it1==nullptr) {
			working_edges[findex].vertices=fbuffer;
			working_edges[findex].first_neighbour=(i/3);

			it2->next_edge=working_edges+findex;
			findex+=1;
			continue;
		}

		if ((it1->vertices).second>fbuffer.second) {
			working_edges[findex].vertices=fbuffer;
			working_edges[findex].first_neighbour=(i/3);

			working_edges[findex].next_edge=it1;

			it2->next_edge=working_edges+findex;
			findex+=1;
			continue;
		}
		if ((it1->vertices).second==fbuffer.second) {
			it1->second_neighbour=(i/3);
			continue;
		}
		std::cout<<"Bad algorithm"<<std::endl;
		exit(1);
		}
	initialized=1;
	nb_edges=findex;
}


void Mesh2d::find_edge(int Simplex_index,std::pair<int,int> edge,int & Ep) {
	Ep=0;
	int index=this->operator()(Simplex_index,Ep);
	while (index==edge.first || index==edge.second) {
		++Ep;
		index=this->operator()(Simplex_index,Ep);
	}
}

int Mesh2d::Adj(int K,int E,int & Ep) {
	if (!initialized) {initialization();}
	std::pair<int,int> fbuffer=wrap_edge(3*K+E);

	Edges * it=register_edges[fbuffer.first];
	while((it->vertices).second<fbuffer.second) {
		it=it->next_edge;
	}
	if (it->first_neighbour==K) {
		if (it->second_neighbour==-1) {
			return -1;
		}
		else {
			find_edge(it->first_neighbour,fbuffer,Ep);
			return it->second_neighbour;
		}
	}
	if (it->second_neighbour==K) {
		find_edge(it->second_neighbour,fbuffer,Ep);
		return it->first_neighbour;
	}
	std::cout<<"Bad algorithm"<<std::endl;
	exit(1);
}


int Mesh2d::IsoK(const R2 * Pk,const double * fhk,R2 * A) {
	if (fhk[0]==0 && fhk[1]==0) {
		A[0]=Pk[0];
		A[1]=Pk[1];
		return 2;
	}
	if (fhk[1]==0 && fhk[2]==0) {
		A[0]=Pk[1];
		A[1]=Pk[2];
		return 2;
	}
	if (fhk[2]==0 && fhk[0]==0) {
		A[0]=Pk[2];
		A[1]=Pk[0];
		return 2;
	}
	if (fhk[0]==0) {
		A[0]=Pk[0];
		return 1;
	}
	if (fhk[1]==0) {
		A[0]=Pk[1];
		return 1;
	}
	if (fhk[2]==0) {
		A[0]=Pk[2];
		return 1;
	}
	if (fhk[0]!=0&&fhk[1]!=0&&fhk[2]!=0) {
		return 0;
	}
	std::cout<<"Bad isovalue"<<std::endl;
	exit(1);
}


double Mesh2d::EikonalK(int n,R2 * A,R* g,R2 P) {
	if (n==2) {
		R2 A0_A1=A[1]-A[0];
		R2 O=proj(P,A0_A1);
		R2 OP=P-O;
		R p=OP.norme();
		R a=A0_A1.norme();
		R lo=-((P-A[0],A[1]-A[0]))/(a*a);
		R d=g[1]-g[0];

		R x=-d*p/(a*std::sqrt(a*a-d*d));
		R l=x+lo;

		if (std::abs(d)<a && 0<=l &&0<=1) {
			return(g[0]-d*lo+d*x+std::sqrt(p*p+x*x));
		}
		R x0=(A[0]-O,A0_A1)/a;
		R x1=(A[1]-O,A0_A1)/a;
		return(std::min(g[0]-d*lo+d*x0+std::sqrt(p*p+x0*x0),
			   		   g[0]-d*lo+d*x1+std::sqrt(p*p+x1*x1)));
	}
	if (n==1) {
		return(g[0]+(P-A[0]).norme());
	}
	std::cout<<"Moment de vérité"<<n<<std::endl;
   assert(n==1 || n==2);
   exit(1);
}

Solve_Eikonale::Solve_Eikonale(const char * m) : Mesh2d(m),u(new double[nv]) {
	initialization();
	queue=new std::tuple<double,int,int>[2*nb_edges];
	FI_queue=0;
 	LI_queue=0;
}

void Solve_Eikonale::queue_push(double v,int Q,int Ep) {
	LI_queue+=1;
	queue[LI_queue]=std::make_tuple(v,Q,Ep);
}

std::tuple<double,int,int> Solve_Eikonale::pop_queue() {
	assert(FI_queue<LI_queue);
	FI_queue+=1;
	return(queue[FI_queue-1]);
}


void Solve_Eikonale::algorithm(R * fhk,R* g) {
	int n=0;
	std::set<int> T0;
	int Ep;

	g[0]=0;
	g[1]=1;

	for (int i=0;i<nv;++i) {
	u[i]=std::numeric_limits<double>::max();
	}

	for (int i=0;i<nt;++i) {
		int a=this->operator()(i,0);
		int b=this->operator()(i,1);
		int c=this->operator()(i,2);
		if (std::min(fhk[a],std::min(fhk[b],fhk[c]))<=0 && std::max(fhk[a],std::max(fhk[b],fhk[c]))>=0) {
			 T0.insert(i);
		}
	}

for (int i=0;i<nv;++i) {
	std::cout<<i<<"\t"<<u[i]<<std::endl;
}

	for (auto it=T0.begin();it!=T0.end();++it) {
		std::cout<<*it<<std::endl;

		R2 isovalues[2];

		int a=this->operator()(*it,0);
		int b=this->operator()(*it,1);
		int c=this->operator()(*it,2);

		R fhk_2[3];
		fhk_2[0]=fhk[a];
		fhk_2[1]=fhk[b];
		fhk_2[2]=fhk[c];

		R2 Pk[3];
		Pk[0]=*t[*it].v[0];
		Pk[1]=*t[*it].v[1];
		Pk[2]=*t[*it].v[2];

		int n=IsoK(Pk,fhk_2,isovalues);
		u[a]=std::min(u[a],EikonalK(n,isovalues,g,Pk[0]));
		u[b]=std::min(u[b],EikonalK(n,isovalues,g,Pk[1]));
		u[c]=std::min(u[c],EikonalK(n,isovalues,g,Pk[2]));
	}
std::cout<<"Afterwards u becomes"<<std::endl;
for (int i=0;i<nv;++i) {
	std::cout<<i<<"\t"<<u[i]<<std::endl;
}
	for (auto it=T0.begin();it!=T0.end();++it) {
		for (int i=0;i<3;++i) {
			int adj=Adj(*it,i,Ep);
			if (adj!=-1 && T0.find(adj)==T0.end()) {
				R2 A[2];
				A[0]=*t[*it].v[(i+1)%3];
				A[1]=*t[*it].v[(i+2)%3];
				R u_[2];
				
				u_[0]=u[this->operator()(*it,(i+1)%3)];
				u_[1]=u[this->operator()(*it,(i+2)%3)];
			
				queue_push(EikonalK(2,A,u_,*t[*it].v[i]),adj,Ep);
			}
		}
	}

std::cout<<"last in ="<<LI_queue<<std::endl;
std::cout<<"first in ="<<FI_queue<<std::endl;


	while (LI_queue!=FI_queue) { //<- equivalent to the queue is non empty 
		std::tuple<double,int,int> first_in=pop_queue();
		if (T0.find(std::get<1>(first_in))==T0.end()) {
			int index=std::get<1>(first_in);
			T0.insert(index);
			u[this->operator()(index,std::get<2>(first_in))]=std::get<0>(first_in);
			for (int i=0;i<3;++i) {
				int adj=Adj(index,i,Ep);
				if (adj!=-1 && T0.find(adj)==T0.end()) {
					R2 A[2];
					A[0]=*t[index].v[(i+1)%3];
					A[1]=*t[index].v[(i+2)%3];
					R u_[2];
				
					u_[0]=u[this->operator()(index,(i+1)%3)];
					u_[1]=u[this->operator()(index,(i+2)%3)];
				
					queue_push(EikonalK(2,A,u_,*t[index].v[i]),adj,Ep);
				}
			}
		}
	}
std::cout<<"last in ="<<LI_queue<<std::endl;
std::cout<<"first in ="<<FI_queue<<std::endl;

for (int i=0;i<nv;++i) {
	std::cout<<"u("<<i<<")=\t"<<u[i]<<std::endl;
}
	sign(fhk);
}

void Solve_Eikonale::sign(R* fhk) {
	for (int i=0;i<nv;++i) {
		if (fhk[i]<0) {
			u[i]=std::abs(u[i]);
		}
		else {
			u[i]=-std::abs(u[i]);
		}
	}
}
