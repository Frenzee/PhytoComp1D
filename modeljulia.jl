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


function RunProgram()
PhytoParams= CSV.read("/home/michael/JuliaCoding/PhytoSpeciesParams.csv")
LakeParams=CSV.read("/home/michael/JuliaCoding/LakeParams.csv")
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
DZ=10.0
S=zeros(N+2)
S[1]=0
S[2:N+1]=[(i-1/2)*DZ for i in 2:N+1]
S[N+2]=1000.0



function phyto1d(du,u,p,t,J)
 for i=2:(MZ+1)
    for j=1:NP
        du[i,j-1]=growthrates[i,j-1]-(J[i,j]-J[i,j-1])/Dzs[j]
    end
 end
end

function calcgrowths!(u,lights,nutrients,temperatures,growthrates,p)
   for i=1:NS
      for j=1:MZ
         growthrates[j,i]=growthrate(lights[j],nutrients[j],temperatures[j],p[i])
      end
   end
end

function growthrate(light::Float64,nutrient::Float64,temperature::Float64,param::PhytoParam)
   out=param.Pmax*(min(light/(light+param.Ik),nutrient/(nutrient+param.Nk)))
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




function calcadvections(ws,u,p)
   uz=u[i,j]
   uz_down= (j == MZ) ? ZERO : u[i,j+1];
   uz_up=(j<=1) ? ZERO : u[i,j-1]);
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

function calcdiffusions(diffs::Array{Float64,1})
   D=diffs[j];
   DUp= (j>=2)? Ith(diffs,j-1) : ZERO
   hdiff=(j!=MZ) ? D*(uz_down-uz)/(Ith(DZs,j)) : ZERO
   hdiff_minus=(j==1) ? ZERO : DUp*(uz-uz_up)/(Ith(DZs,j))
end

function calcJs(J::Array{Float64,2})
 J[1]=0
 J[MZ+1]=0
 for i=2:NP
    for i=2:MZ
      if (ws<0 & i==MZ)
         J[i]=ws*((uz_down+uz)*HALF)
      elseif ((ws>0) & (i==2))
         J[i]=ws*((uz_down+uz)*HALF)
      elseif (ws<0)
         J[i]=hadv[i-1]-hdiff[i-1]
      elseif (ws>0)
      end
   end
end
end
