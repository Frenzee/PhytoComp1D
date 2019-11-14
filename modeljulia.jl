using Pkg
using DifferentialEquations
using Plots
using CSV
using Flux
using DiffEqFlux
const MZ=100
const NP=10
const NN=0
const NZ=0
const NS=NP+NN
const NEQ=NS*MZ
const DZ=0.1
const SIXTH=1.0/6.0
const NOPARAMS=16

struct PhytoParam
   Pmax::Float64
   Ik::Float64
   ks::Float64
   l::Float64
   Nk::Float64
   a::Float64
   A::Float64
   b::Float64
   B::Float64
   ws::Float64
   u::Array{Float64,1}
end

struct Nutrient
   u::Array{Float64,1}
   ws::Float64

end

struct Virus
   u::Array{Float,1}
   p1::Float64
end

struct Zooplankton

end

HOMEPATH="/home/michael/JuliaCoding/"
function RunProgram()
PhytoParams= CSV.read("HOMEPATH/PhytoSpeciesParams.csv")
LakeParams=CSV.read("HOMEPATH/LakeParams.csv")
lights=zeros(Float64,MZ)
calclights!(u,Iin,Kbg,ks,I,Dzs)
end

#==
1. Calculate light field
2. Calculate growth rates, diffusions and advections
4. Populate differential equation matrix
5. Solve for one timestep
6. GOTO 1.
==#
N=100
DZ=10
NS=10
S=zeros(Float64,N+2)
S[1]=0
S[2:N+1]=[(i-1/2)*DZ for i in 2:N+1]
S[N+2]=1000.0

Js=zeros(Float64,N+1,NS)
for i in 2:N
   for j in 1:NS
      Js[i,j]= (i==2)? v*(w[1,j]+w[2,j])/2.0-D*(w[2,j]-w[1,j])/DZ[i] : v*1.0/6.0(-w[i-1,j]+5*w[i,j]+2*w[i+1,j])
end
end


function calcgrowths!(u,lights,nutrients,temperatures,growthrates,p)
   for i=1:NS
      for j=1:MZ
         growthrates[j,i]=growthrate(lights[j],nutrients[j],temperatures[j],p)
      end
   end
end

function growthrate(light::Float64,nutrient::Float64,temperature::Float64,param::PhytoParam)
   out=param.Pmax*(min(light/(light+param.Ik),nutrient/(nutrient+param.Nk)))
end

function temperaturegrowth(A::Float64,a::Float64,B::Float64,b::Float64,temperature::Float64)
   out=A*exp(a*temperature)-B*exp(b*temperature)
end

function calcweightedbiomass(u::Array{Float64,2},ks::Array{Float64,1})
   mz=length(u[1,:])
   np=length(ks)
   out=zeros(Float64,mz)
   for i=1:mz
      for j=1:np
         out[i]+=u[j,i]*ks[j]
      end
   end
   return(out)
end

function calclights!(u::Array{Float64,2},Iin::Float64,Kbg::Float64,ks::Array{Float64,1},Is::Array{Float64,1},Dzs::Array{Float64,1})
   uSum=calcweightedbiomass(u,ks)
   ws=Array{Float64,1}(undef,length(u[:,1]))
   w0=(3.0*uSum[1]-uSum[2])/2.0
   ws[1]=1.0/2.0*(w0+uSum[1]+Kbg)*Dzs[1]
   Is[1]=Iin*exp(-ws[1])
   i=2
   while i<=length(ws)
      ws[i]=ws[i-1]+(1.0/2.0*(uSum[i-1]+uSum[i])+Kbg)*Dzs[i]
      Is[i]=Iin*exp(-ws[i])
      i+=1
   end
   Is
end

function calcadvections(ws,u,p,i,j)
   uz=u[i,j]
   uz_down= (j == MZ) ? ZERO : u[i,j+1];
   uz_up=(j<=1) ? ZERO : u[i,j-1];
   ws=PhytoParams.ws[i]
   if ws<0
      uz_downtwo= (j >MZ-2) ? ZERO : u[i,j+2];
      hadv=p.ws*(SIXTH)*(TWO*uz + FIVE*uz_down - uz_downtwo)
      hadv_minus=p.ws*(SIXTH)*(TWO*uz_up + FIVE*uz -uz_down)
   else
      uz_uptwo=(j<=2) ? ZERO : u[i,j-2];
      hadv=p.ws*SIXTH*(-uz_up+FIVE*uz+TWO*uz_down)
      hadv_minus= p.ws*SIXTH*(-uz_uptwo+FIVE*uz_up+TWO*uz)
   end
end

function calcdiffusions(diffs::Array{Float64,1},u,i,j)
   D=diffs[j];
   DUp= (j>=2)? Ith(diffs,j-1) : ZERO
   hdiff=(j!=MZ) ? D*(uz_down-uz)/(Ith(DZs,j)) : ZERO
   hdiff_minus=(j==1) ? ZERO : DUp*(uz-uz_up)/(Ith(DZs,j))
end
