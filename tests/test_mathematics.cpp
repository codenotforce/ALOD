#include "alod/problems.hpp"
#include "alod/meshes.hpp"
#include "alod/mesh_state.hpp"
#include "alod/afem.hpp"
#include "helmholtz/boundary.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>

using namespace lod2d;using namespace lod2d::helmholtz;
void require(bool ok,const char* text){if(!ok)throw std::runtime_error(text);}
template<class F> void rejects(F f){bool caught=false;try{f();}catch(const std::exception&){caught=true;}require(caught,"invalid input was accepted");}

void manufactured(std::string id) {
    auto p=alod::make_problem(id);
    std::vector<Point2> points=id=="E1"?std::vector<Point2>{{.71,.43},{.5,.51},{.83,.64}}:
        std::vector<Point2>{{-.48,.52},{-.1,.08},{.2,.17},{-.4,-.3}};
    double gradient_error=0,laplacian_error=0;
    for(auto x:points){
        require(std::abs(p.source(x)+p.exact_laplacian(x)+256.0*p.exact(x))<1e-10*(1+std::abs(p.source(x))),"manufactured PDE identity");
        double h=1e-5;Eigen::Vector2cd g;Complex lap=0;
        for(int d=0;d<2;++d){Point2 a=x,b=x;a[d]+=h;b[d]-=h;
            g[d]=(p.exact(a)-p.exact(b))/(2*h);
            lap+=(p.exact_gradient(a)[d]-p.exact_gradient(b)[d])/(2*h);
        }
        gradient_error=std::max(gradient_error,(g-p.exact_gradient(x)).norm()/(1+p.exact_gradient(x).norm()));
        laplacian_error=std::max(laplacian_error,std::abs(lap-p.exact_laplacian(x))/(1+std::abs(p.exact_laplacian(x))));
    }
    require(gradient_error<3e-7 && laplacian_error<3e-7,"independent manufactured derivatives");
    double boundary_error=0;
    for(const auto& edge:p.initial_mesh.boundary_edges){
        Point2 a=p.initial_mesh.nodes[edge.nodes[0]],b=p.initial_mesh.nodes[edge.nodes[1]],v=b-a;
        Point2 n(v.y(),-v.x());n.normalize();
        for(const auto& t:p.initial_mesh.elems) {
            bool has_a=std::find(t.begin(),t.end(),edge.nodes[0])!=t.end();
            bool has_b=std::find(t.begin(),t.end(),edge.nodes[1])!=t.end();
            if(has_a && has_b)for(int i:t)if(i!=edge.nodes[0] && i!=edge.nodes[1])
                if(n.dot(p.initial_mesh.nodes[i]-(a+b)/2)>0)n=-n;
        }
        for(double t:{.15,.37,.63,.85}){Point2 x=(1-t)*a+t*b;Complex r=p.exact(x);
            if(edge.tag!=BoundaryTag::Dirichlet)r=n.cast<Complex>().dot(p.exact_gradient(x))-(edge.tag==BoundaryTag::Robin?Complex(0,16)*p.exact(x):Complex(0));
            boundary_error=std::max(boundary_error,std::abs(r));
        }
    }
    require(boundary_error<1e-11,"manufactured mixed boundary residual");
    auto mesh=refine_mesh_nvb(p.initial_mesh,2).mesh;
    auto energy=[&](QuadraturePolicy q){return integrate_scalar_function(mesh,[&](const Point2& x){return p.exact_gradient(x).squaredNorm()+256*std::norm(p.exact(x));},q,p.quadrature_context);};
    double low=energy({4,4,6,1}),mid=energy({8,12,16,4});
    double high=energy(alod::paper_quadrature(id)),reference=energy({18,26,32,10});
    std::cout<<id<<" derivative="<<gradient_error<<","<<laplacian_error<<" boundary="<<boundary_error
        <<" quadrature="<<low<<","<<mid<<","<<high<<","<<reference<<'\n';
    require(std::abs(high-reference)<1e-6*reference,"paper quadrature accuracy");
    require(std::abs(high-reference)<std::abs(low-reference),"quadrature convergence from low order");
    require(std::abs(mid-reference)<std::abs(low-reference),"intermediate quadrature convergence");
}

void meshes(std::string id){
    auto p=alod::make_problem(id);alod::MeshState a(p.initial_mesh),b(p.initial_mesh);
    auto boundary=[&](BoundaryTag tag){return boundary_measure(p.initial_mesh,tag);};
    std::set<std::uint64_t> known;for(auto x:a.ancestry)known.insert(x.id);
    for(int step=0;step<6;++step){
        auto old=a.mesh;auto ids=a.elements;
        std::vector<int> marks{static_cast<int>((step*3)%old.elems.size())};
        auto out=bisect_newest_vertex(old,marks);
        Eigen::VectorXd linear(old.nodes.size());for(int i=0;i<linear.size();++i)linear[i]=2+old.nodes[i].x()-3*old.nodes[i].y();
        Eigen::VectorXd injected=out.P_node*linear;
        for(int i=0;i<injected.size();++i)require(std::abs(injected[i]-(2+out.mesh.nodes[i].x()-3*out.mesh.nodes[i].y()))<1e-13,"NVB injection");
        a.refine(marks);b.refine(marks);validate_boundary_tags(a.mesh);
        require(alod::mesh_fingerprint(a.mesh)==alod::mesh_fingerprint(b.mesh),"deterministic mesh fingerprint");
        require(a.elements.size()==b.elements.size(),"stable identity size");
        for(std::size_t i=0;i<a.elements.size();++i){auto x=a.elements[i];
            require(x.id==b.elements[i].id && x.parent==b.elements[i].parent,"stable identities and parents");
            if(!known.contains(x.id)){require(known.contains(x.parent),"new child has no prior parent");known.insert(x.id);}
        }
        for(auto tag:{BoundaryTag::Dirichlet,BoundaryTag::Neumann,BoundaryTag::Robin})
            require(std::abs(boundary_measure(a.mesh,tag)-boundary(tag))<1e-13,"NVB boundary measure");
        auto areas=compute_area(a.mesh);double area=std::accumulate(areas.begin(),areas.end(),0.0);
        require(*std::min_element(areas.begin(),areas.end())>0 && std::abs(area-(id=="E1"?1:3))<1e-13,"NVB area conservation");
        const auto [edges,is_boundary]=compute_edges(a.mesh);
        for(auto e:edges){Point2 x=a.mesh.nodes[e[0]],d=a.mesh.nodes[e[1]]-x;
            for(int i=0;i<static_cast<int>(a.mesh.nodes.size());++i)if(i!=e[0]&&i!=e[1]){
                Point2 v=a.mesh.nodes[i]-x;double t=v.dot(d)/d.squaredNorm();
                require(!(t>1e-12&&t<1-1e-12&&std::abs(v.x()*d.y()-v.y()*d.x())<1e-13),"hanging node");
            }
        }
    }
    auto fingerprint=alod::mesh_fingerprint(a.mesh);auto version=a.version;a.refine({});
    require(fingerprint==alod::mesh_fingerprint(a.mesh)&&version==a.version,"empty refinement must be identity");
    rejects([&]{a.refine({-1});});
}

void residuals(){
    auto mesh=make_helmholtz_unit_square_mesh();auto op=assemble_helmholtz_operators(mesh,2);
    ComplexVector u=ComplexVector::Ones(mesh.nodes.size());
    ComplexFunction f=[](const Point2&){return Complex(-4,0);};
    auto b=assemble_helmholtz_load(mesh,f);
    auto r=adaptive::diagnostics::estimate_conforming_p1_residual(mesh,op,u,b,f);
    double total=std::accumulate(r.element_squared.begin(),r.element_squared.end(),0.0);
    require(std::abs(total-16)<1e-11,"Robin edge integration must equal k^2 times squared edge lengths");
    require(r.algebraic_relative_difference<1e-12,"assembled residual matches independent strong residual");
    for(double x:r.interior_jump_squared)require(x==0,"constant field has no interior jump");
    require(adaptive::mark_doerfler({1,1,1,1},.5)==std::vector<int>({0,1}),"deterministic Doerfler tie break");
    require(adaptive::mark_doerfler({1,4,2,3},.5)==std::vector<int>({1,3}),"Doerfler bulk uses squared indicators");
    rejects([]{adaptive::mark_doerfler({1,-1},.5);});
    rejects([]{adaptive::mark_doerfler({1,2},0);});
    rejects([]{adaptive::mark_doerfler({1,9},.5,{1,0});});
}
int main(){try{for(std::string id:{"E1","E2"}){manufactured(id);meshes(id);}residuals();
    rejects([]{alod::make_problem("E3");});std::cout<<"mathematical and mesh checks passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
