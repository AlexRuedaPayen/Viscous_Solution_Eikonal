#include "EF2d-base.hpp"
#include <iostream>
#include <string>
#include <ctime>
#include <cmath>
#include <algorithm>

double distance(R2 input) {
	return(input.y-0.5);
}

double hole_in_mesh(R2 input) {
	return(std::sqrt(input.x*input.x+input.y*input.y)-1);
}

double double_elliptic_hole(double * x) {
	int a=1/1.33;
	int b=1/0.61;
	double c=std::sqrt(a*a*x[0]*x[0]+x[1]*x[1]);
	double d=std::sqrt(x[0]*x[0]+b*b*x[1]*x[1]);
	return(c*d);
}

void help() {
	std::cout<<"\n";
	std::cout<<"Hello this interface is dedicated to help"<<std::endl;
	std::cout<<"It is still in developpement"<<std::endl;
	std::cout<<"\n";
	std::cout<<"OPTIONS :\n";
	std::cout<<"\t-solve [meshname]\n";
	std::cout<<"\t\tReads [filename] and solve the Eikonal problem\n";
	std::cout<<"\n";
	std::cout<<"\t-debug [meshname]\n";
	std::cout<<"\t\tReads [filename] and carry out tests to check the function developed for the project\n";
	std::cout<<"\n";
	std::cout<<"Bye, Thank you for using this algorithm"<<std::endl<<std::endl;
	exit(0);
}

void solve(int &CPUmemspace,const char * filename) {
	Solve_Eikonale Th(filename);
	Th.initialization();
	double fhk_[Th.nv];
	if ((std::string)filename=="./Mesh/distance.msh") {
		for (int i=0;i<Th.nv;++i) {
			fhk_[i]=distance(Th.v[i]);
			std::cout<<fhk_[i]<<std::endl;
		}
	}
	if((std::string)filename=="./Mesh/Thg.msh") {
		for (int i=0;i<Th.nv;++i) {
			fhk_[i]=distance(Th.v[i]);
			std::cout<<fhk_[i]<<std::endl;
		}
	}
	R g[2];
	g[0]=0;
	g[1]=0;
	/*
	Th.algorithm(fhk_,g);
	CPUmemspace=std::max(CPUmemspace,(int)sizeof(Th));
	*/
}

void debug(int &CPUmemspace,const char * filename) {
	Mesh2d Thp(filename);
	std::cout<<"nv = "<<Thp.nv<<std::endl;
	std::cout<<"nt = "<<Thp.nt<<std::endl;
	for (int i=0;i<Thp.nt;++i) {
		std::cout<<Thp(i,0)<<"\t"<<Thp(i,1)<<"\t"<<Thp(i,2)<<"\t"<<std::endl;
	}
	Thp.initialization();
	for (int i=0;i<(3*Thp.nt);++i) {
		if (Thp.working_edges+i!=nullptr) {
			std::cout<<"("<<Thp.working_edges[i].vertices.first<<","<<Thp.working_edges[i].vertices.second<<")";
			std::cout<<"\t"<<Thp.working_edges[i].first_neighbour<<"\t"<<Thp.working_edges[i].second_neighbour<<std::endl;
		}
	}
	int Ep;
	
	for (int triangle_test=0;triangle_test<Thp.nt;++triangle_test) {
		int test_Adj=Thp.Adj(triangle_test,0,Ep);
		std::cout<<test_Adj<<std::endl;
		if (test_Adj!=-1) {
			std::cout<<Thp(test_Adj,0)<<"\t"<<Thp(test_Adj,1)<<"\t"<<Thp(test_Adj,2)<<"\n";
			std::cout<<Thp(triangle_test,0)<<"\t"<<Thp(triangle_test,1)<<"\t"<<Thp(triangle_test,2)<<"\n";
			std::cout<<Ep<<"\n";
		}
		test_Adj=Thp.Adj(triangle_test,1,Ep);
		std::cout<<test_Adj<<std::endl;
		if (test_Adj!=-1) {
			std::cout<<Thp(test_Adj,0)<<"\t"<<Thp(test_Adj,1)<<"\t"<<Thp(test_Adj,2)<<"\n";
			std::cout<<Thp(triangle_test,0)<<"\t"<<Thp(triangle_test,1)<<"\t"<<Thp(triangle_test,2)<<"\n";
			std::cout<<Ep<<"\n";
		}
		test_Adj=Thp.Adj(triangle_test,2,Ep);
		std::cout<<test_Adj<<std::endl;
		if (test_Adj!=-1) {
			std::cout<<Thp(test_Adj,0)<<"\t"<<Thp(test_Adj,1)<<"\t"<<Thp(test_Adj,2)<<"\n";
			std::cout<<Thp(triangle_test,0)<<"\t"<<Thp(triangle_test,1)<<"\t"<<Thp(triangle_test,2)<<"\n";
			std::cout<<Ep<<"\n";
		}
	}
	CPUmemspace=sizeof(Thp);

	R2 Triangle[3];
	Triangle[0]=R2(1.5,2.3);
	Triangle[1]=R2(0.8,-3);
	Triangle[2]=R2(-0.2,0.4);
	double fhk[3];
	fhk[0]=-1;
	fhk[1]=3.14;
	fhk[2]=1.1618;
	R2 p[2];
	
	
	std::cout<<"IsoK\t"<<Thp.IsoK(Triangle,fhk,p)<<std::endl;

	R2 proj_test=proj(R2(1,-1),R2(1,5));
	std::cout<<proj_test.x<<"\t"<<proj_test.y<<"\n";

	R g[2];
	g[0]=0;
	g[1]=0;
	double fhk_[Thp.nv];

	if ((std::string)filename=="./Mesh/distance.msh") {
		for (int i=0;i<Thp.nv;++i) {
			fhk_[i]=distance(Thp.v[i]);
			std::cout<<fhk_[i]<<std::endl;
		}
	}
	if((std::string)filename=="./Mesh/Thg.msh") {
		for (int i=0;i<Thp.nv;++i) {
			fhk_[i]=distance(Thp.v[i]);
			std::cout<<fhk_[i]<<std::endl;
		}
	}
/*
	Solve_Eikonale Eik_Thp(filename);
	Eik_Thp.algorithm(fhk_,g);
	CPUmemspace=std::max(CPUmemspace,(int)sizeof(Eik_Thp));
*/
}


int main(int argc,const char ** argv) {
	clock_t start,end;
	start=clock();
	int CPUmemspace=0;
	if (argc>1) {
		std::string current_option_name=argv[1];
		if (current_option_name=="-debug") {
			assert(argc>2);
			debug(CPUmemspace,argv[2]);
		}
		if (current_option_name=="-help") {
			help();
		}
		if (current_option_name=="-solve") {
			assert(argc>2);
			solve(CPUmemspace,argv[2]);
		}
	}
	end=clock();
	double time_in_seconds=static_cast<double>(end-start) /CLOCKS_PER_SEC;
	std::cout<<"\n"<<std::endl;
	std::cout<<"Algorithms written by Alexandre Rueda Payen and Paul Hervé from Frederic Hecht's preliminary work (2020)"<<std::endl;
	std::cout<<"\n"<<std::endl;
	std::cout<<"Time taken :\t"<<time_in_seconds<<" seconds"<<std::endl;
	std::cout<<"Maximal CPU memory used :\t"<<CPUmemspace<<" bytes"<<std::endl;
	std::cout<<"\n"<<std::endl;
	std::cout<<"Thank you for running this algorithm ^_^"<<std::endl;
	std::cout<<"\n"<<std::endl;
	return 0;
}