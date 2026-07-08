/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-2002 Vincent Favre-Nicolin vincefn@users.sourceforge.net
        2000-2001 University of Geneva (Switzerland)

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*
*  source file ObjCryst++ Crystal class
*
*/
#include "ObjCryst/ObjCryst/UnitCell.h"
#include "ObjCryst/Quirks/VFNStreamFormat.h"

namespace ObjCryst
{
const RefParType *gpRefParTypeUnitCell=0;
const RefParType *gpRefParTypeUnitCellLength=0;
const RefParType *gpRefParTypeUnitCellAngle=0;
long NiftyStaticGlobalObjectsInitializer_UnitCell::mCount=0;

UnitCell::UnitCell():
mCellDim(6),
mBMatrix(3,3),mOrthMatrix(3,3),mOrthMatrixInvert(3,3)
{
   VFN_DEBUG_MESSAGE("UnitCell::UnitCell()",10)
   this->InitOptions();
   this->Init(10,11,12,M_PI/2+.1,M_PI/2+.2,M_PI/2+.3,"P 1","");
   mClockMaster.AddChild(mClockLatticePar);
   mClockMaster.AddChild(this->GetSpaceGroup().GetClockSpaceGroup());
}

UnitCell::UnitCell(const REAL a, const REAL b, const REAL c, const string &SpaceGroupId):
mCellDim(6),
mBMatrix(3,3),mOrthMatrix(3,3),mOrthMatrixInvert(3,3)
{
   VFN_DEBUG_MESSAGE("UnitCell::UnitCell(a,b,c,Sg)",10)
   this->Init(a,b,c,M_PI/2,M_PI/2,M_PI/2,SpaceGroupId,"");
   this->InitOptions();
   mClockMaster.AddChild(mClockLatticePar);
   mClockMaster.AddChild(this->GetSpaceGroup().GetClockSpaceGroup());
}

UnitCell::UnitCell(const REAL a, const REAL b, const REAL c, const REAL alpha,
              const REAL beta, const REAL gamma,const string &SpaceGroupId):
mCellDim(6),
mBMatrix(3,3),mOrthMatrix(3,3),mOrthMatrixInvert(3,3)
{
   VFN_DEBUG_MESSAGE("UnitCell::UnitCell(a,b,c,alpha,beta,gamma,Sg)",10)
   this->Init(a,b,c,alpha,beta,gamma,SpaceGroupId,"");
   this->InitOptions();
   mClockMaster.AddChild(mClockLatticePar);
   mClockMaster.AddChild(this->GetSpaceGroup().GetClockSpaceGroup());
}

UnitCell::UnitCell(const UnitCell &old):
mCellDim(old.mCellDim),mSpaceGroup(old.GetSpaceGroup()),
mBMatrix(3,3),mOrthMatrix(3,3),mOrthMatrixInvert(3,3)
{
   VFN_DEBUG_MESSAGE("UnitCell::UnitCell(&oldUnitCell)",10)
   this ->InitRefParList();
   mConstrainLatticeToSpaceGroup.SetChoice(old.mConstrainLatticeToSpaceGroup.GetChoice());
   this->InitMatrices();
   this->UpdateLatticePar();
   mClockMaster.AddChild(mClockLatticePar);
   mClockMaster.AddChild(this->GetSpaceGroup().GetClockSpaceGroup());
}

UnitCell::~UnitCell()
{
   VFN_DEBUG_ENTRY("UnitCell::~UnitCell()",5)
   VFN_DEBUG_EXIT("UnitCell::~UnitCell()",5)
}

const string& UnitCell::GetClassName() const
{
   const static string className="UnitCell";
   return className;
}

CrystVector_REAL UnitCell::GetLatticePar() const
{
   VFN_DEBUG_MESSAGE("UnitCell::GetLatticePar()",0)

   if(  (mClockLatticeParUpdate>mClockLatticePar)
      &&(mClockLatticeParUpdate>mSpaceGroup.GetClockSpaceGroup())) return mCellDim;
   else
   {
      //:NOTE: cannot use this->UpdateLatticePar() because it is not a const member function
      int num = mSpaceGroup.GetSpaceGroupNumber();
      CrystVector_REAL cellDim=mCellDim;
      if((num <=2)||(mConstrainLatticeToSpaceGroup.GetChoice()!=0)) return cellDim;
      if((num <=15) && (0==mSpaceGroup.GetUniqueAxis()))
      {
         cellDim(4)=M_PI/2.;
         cellDim(5)=M_PI/2.;
         return cellDim;
      }
      if((num <=15) && (1==mSpaceGroup.GetUniqueAxis()))
      {
         cellDim(3)=M_PI/2.;
         cellDim(5)=M_PI/2.;
         return cellDim;
      }
      if((num <=15) && (2==mSpaceGroup.GetUniqueAxis()))
      {
         cellDim(3)=M_PI/2.;
         cellDim(4)=M_PI/2.;
         return cellDim;
      }
      if(num <=74)
      {
         cellDim(3)=M_PI/2.;
         cellDim(4)=M_PI/2.;
         cellDim(5)=M_PI/2.;
         return cellDim;
      }
      if(num <= 142)
      {
         cellDim(3)=M_PI/2.;
         cellDim(4)=M_PI/2.;
         cellDim(5)=M_PI/2.;
         cellDim(1) = mCellDim(0) ;
         return cellDim;
      }
      if(mSpaceGroup.GetExtension()=='R')
      {
         cellDim(4) = mCellDim(3);
         cellDim(5) = mCellDim(3);
         cellDim(1) = mCellDim(0);
         cellDim(2) = mCellDim(0);
         return cellDim;
      }
      if(num <= 194) // || (mSpaceGroup.GetExtension()=='H') //Hexagonal
      {//Hexagonal axes, for hexagonal and non-rhomboedral trigonal cells
         cellDim(3) = M_PI/2.;
         cellDim(4) = M_PI/2.;
         cellDim(5) = M_PI*2./3.;
         cellDim(1) = mCellDim(0) ;
         return cellDim;
      }
      cellDim(3)=M_PI/2.;
      cellDim(4)=M_PI/2.;
      cellDim(5)=M_PI/2.;
      cellDim(1) = mCellDim(0) ;
      cellDim(2) = mCellDim(0) ;
      return cellDim;
   }
}

REAL UnitCell::GetLatticePar(int whichPar)const
{
   VFN_DEBUG_MESSAGE("UnitCell::LatticePar(i)",0)
   if( (whichPar<0) || (whichPar>5))
      throw ObjCrystException("UnitCell::LatticePar(int) :trying to access parameter>5!");

   if(  (mClockLatticeParUpdate>mClockLatticePar)
      &&(mClockLatticeParUpdate>mSpaceGroup.GetClockSpaceGroup())) return mCellDim(whichPar);
   else
   {
      const int num = mSpaceGroup.GetSpaceGroupNumber();

      static CrystVector_REAL cellDim;
      cellDim=mCellDim;
      if((num <=2)||(mConstrainLatticeToSpaceGroup.GetChoice()!=0))
         return cellDim(whichPar);
      if((num <=15) && (0==mSpaceGroup.GetUniqueAxis()))
      {
         cellDim(4)=M_PI/2.;
         cellDim(5)=M_PI/2.;
         return cellDim(whichPar);
      }
      if((num <=15) && (1==mSpaceGroup.GetUniqueAxis()))
      {
         cellDim(3)=M_PI/2.;
         cellDim(5)=M_PI/2.;
         return cellDim(whichPar);
      }
      if((num <=15) && (2==mSpaceGroup.GetUniqueAxis()))
      {
         cellDim(3)=M_PI/2.;
         cellDim(4)=M_PI/2.;
         return cellDim(whichPar);
      }

      if(num <=74)
      {
         cellDim(3)=M_PI/2.;
         cellDim(4)=M_PI/2.;
         cellDim(5)=M_PI/2.;
         return cellDim(whichPar);
      }
      if(num <= 142)
      {
         cellDim(3)=M_PI/2.;
         cellDim(4)=M_PI/2.;
         cellDim(5)=M_PI/2.;
         cellDim(1) = mCellDim(0) ;
         return cellDim(whichPar);
      }
      if(mSpaceGroup.GetExtension()=='R')
      {
         cellDim(4) = mCellDim(3);
         cellDim(5) = mCellDim(3);
         cellDim(1) = mCellDim(0);
         cellDim(2) = mCellDim(0);
         return cellDim(whichPar);
      }
      if(num <= 194) // ||(mSpaceGroup.GetExtension()=='H')
      {//Hexagonal axes, for hexagonal and non-rhomboedral trigonal cells
         cellDim(3) = M_PI/2.;
         cellDim(4) = M_PI/2.;
         cellDim(5) = M_PI*2./3.;
         cellDim(1) = mCellDim(0) ;
         return cellDim(whichPar);
      }
      cellDim(3)=M_PI/2.;
      cellDim(4)=M_PI/2.;
      cellDim(5)=M_PI/2.;
      cellDim(1) = mCellDim(0) ;
      cellDim(2) = mCellDim(0) ;
      return cellDim(whichPar);
   }
}

const CrystMatrix_REAL& UnitCell::GetBMatrix()const
{
   VFN_DEBUG_MESSAGE("UnitCell::GetBMatrix()",0)
   this->InitMatrices();
   return mBMatrix;
}

const CrystMatrix_REAL& UnitCell::GetOrthMatrix() const
{
   VFN_DEBUG_MESSAGE("UnitCell::GetOrthMatrix()",0)
   this->InitMatrices();
   return mOrthMatrix;
}

CrystVector_REAL UnitCell::GetOrthonormalCoords(const REAL x,
                                                const REAL y,
                                                const REAL z) const
{
   this->InitMatrices();
   CrystVector_REAL coords(3);
   coords(0)=mOrthMatrix(0,0)*x+mOrthMatrix(0,1)*y+mOrthMatrix(0,2)*z;
   coords(1)=mOrthMatrix(1,0)*x+mOrthMatrix(1,1)*y+mOrthMatrix(1,2)*z;
   coords(2)=mOrthMatrix(2,0)*x+mOrthMatrix(2,1)*y+mOrthMatrix(2,2)*z;
   return coords;
}

void UnitCell::FractionalToOrthonormalCoords(REAL &x,REAL &y,REAL &z) const
{
   this->InitMatrices();
   const REAL oldx=x;
   const REAL oldy=y;
   x=mOrthMatrix(0,0)*oldx+mOrthMatrix(0,1)*oldy+mOrthMatrix(0,2)*z;
   y=mOrthMatrix(1,0)*oldx+mOrthMatrix(1,1)*oldy+mOrthMatrix(1,2)*z;
   z=mOrthMatrix(2,0)*oldx+mOrthMatrix(2,1)*oldy+mOrthMatrix(2,2)*z;
}

void UnitCell::OrthonormalToFractionalCoords(REAL &x,REAL &y,REAL &z) const
{
   //cout << x << " " << y << " " << z <<endl;
   //cout << endl << mOrthMatrixInvert <<endl;
   this->InitMatrices();
   const REAL oldx=x;
   const REAL oldy=y;
   x=mOrthMatrixInvert(0,0)*oldx+mOrthMatrixInvert(0,1)*oldy+mOrthMatrixInvert(0,2)*z;
   y=mOrthMatrixInvert(1,0)*oldx+mOrthMatrixInvert(1,1)*oldy+mOrthMatrixInvert(1,2)*z;
   z=mOrthMatrixInvert(2,0)*oldx+mOrthMatrixInvert(2,1)*oldy+mOrthMatrixInvert(2,2)*z;
   //cout << x << " " << y << " " << z <<endl;
}

void UnitCell::MillerToOrthonormalCoords(REAL &x,REAL &y,REAL &z) const
{
   this->InitMatrices();
   const REAL oldx=x;
   const REAL oldy=y;
   x=mBMatrix(0,0)*oldx+mBMatrix(0,1)*oldy+mBMatrix(0,2)*z;
   y=mBMatrix(1,0)*oldx+mBMatrix(1,1)*oldy+mBMatrix(1,2)*z;
   z=mBMatrix(2,0)*oldx+mBMatrix(2,1)*oldy+mBMatrix(2,2)*z;
}

void UnitCell::OrthonormalToMillerCoords(REAL &x,REAL &y,REAL &z) const
{
   this->InitMatrices();
   const REAL oldx=x;
   const REAL oldy=y;
   x=mBMatrixInvert(0,0)*oldx+mBMatrixInvert(0,1)*oldy+mBMatrixInvert(0,2)*z;
   y=mBMatrixInvert(1,0)*oldx+mBMatrixInvert(1,1)*oldy+mBMatrixInvert(1,2)*z;
   z=mBMatrixInvert(2,0)*oldx+mBMatrixInvert(2,1)*oldy+mBMatrixInvert(2,2)*z;
}

void UnitCell::Print(ostream &os)const
{
   VFN_DEBUG_MESSAGE("UnitCell::Print()",5)
   this->InitMatrices();
   int width =8 ;
   os << "UnitCell : " << mName <<"("<<this->GetSpaceGroup().GetName()<<")"<< endl;
   os.width(width);
   os   << "    Cell dimensions : "
        << FormatFloat(this->GetLatticePar(0),8,5) << "  "
        << FormatFloat(this->GetLatticePar(1),8,5) << "  "
        << FormatFloat(this->GetLatticePar(2),8,5) << "  "
        << FormatFloat(this->GetLatticePar(3)*RAD2DEG,8,5) << "  "
        << FormatFloat(this->GetLatticePar(4)*RAD2DEG,8,5) << "  "
        << FormatFloat(this->GetLatticePar(5)*RAD2DEG,8,5) << endl ;
}

const SpaceGroup & UnitCell::GetSpaceGroup() const {return mSpaceGroup;}
SpaceGroup & UnitCell::GetSpaceGroup()  {return mSpaceGroup;}

void UnitCell::ChangeSpaceGroup(const string &spgId)
{
  this->GetSpaceGroup().ChangeSpaceGroup(spgId);
  this->InitRefParList();
  this->UpdateLatticePar();
}


const RefinableObjClock& UnitCell::GetClockLatticePar()const {return mClockLatticePar;}
const RefinableObjClock& UnitCell::GetClockMetricMatrix()const {return mClockMetricMatrix;}

REAL UnitCell::GetVolume()const
{
   const REAL a=this->GetLatticePar(0);
   const REAL b=this->GetLatticePar(1);
   const REAL c=this->GetLatticePar(2);
   const REAL alpha=this->GetLatticePar(3);
   const REAL beta=this->GetLatticePar(4);
   const REAL gamma=this->GetLatticePar(5);

   return a*b*c*sqrt(1-cos(alpha)*cos(alpha)-cos(beta)*cos(beta)-cos(gamma)*cos(gamma)
            +2*cos(alpha)*cos(beta)*cos(gamma));
}

CrystMatrix_REAL UnitCell::GetReciprocalMetricTensorDeriv(const RefinablePar &par) const
{
   CrystMatrix_REAL dGstar(3,3);
   dGstar=0;
   // Which stored lattice parameter is this ?
   int idx=-1;
   for(int i=0;i<6;i++) if(par.GetPointer()==(mCellDim.data()+i)) {idx=i;break;}
   if(idx<0) return dGstar;

   // Constraint mapping: which *effective* lattice parameters (as returned by
   // GetLatticePar()) change when this stored parameter changes ? Must follow
   // exactly the constraints applied in UnitCell::GetLatticePar().
   bool vAffected[6]={false,false,false,false,false,false};
   vAffected[idx]=true;
   const int num=mSpaceGroup.GetSpaceGroupNumber();
   if((num>2)&&(mConstrainLatticeToSpaceGroup.GetChoice()==0))
   {
      if(num<=15)
      {// Monoclinic: two angles are fixed to pi/2
         const int ax=mSpaceGroup.GetUniqueAxis();
         if((idx>=3)&&(idx!=(ax+3))) vAffected[idx]=false;
      }
      else if(num<=74)
      {// Orthorhombic
         if(idx>=3) vAffected[idx]=false;
      }
      else if(num<=142)
      {// Tetragonal: b=a, c free
         if(idx==0) vAffected[1]=true;
         if((idx==1)||(idx>=3)) vAffected[idx]=false;
      }
      else if(mSpaceGroup.GetExtension()=='R')
      {// Rhombohedral: b=c=a, beta=gamma=alpha
         if(idx==0) {vAffected[1]=true;vAffected[2]=true;}
         else if(idx==3) {vAffected[4]=true;vAffected[5]=true;}
         else vAffected[idx]=false;
      }
      else if(num<=194)
      {// Trigonal/hexagonal (hexagonal axes): b=a, c free
         if(idx==0) vAffected[1]=true;
         if((idx==1)||(idx>=3)) vAffected[idx]=false;
      }
      else
      {// Cubic: b=c=a
         if(idx==0) {vAffected[1]=true;vAffected[2]=true;}
         if(idx>=1) vAffected[idx]=false;
      }
   }
   if(!(vAffected[0]||vAffected[1]||vAffected[2]||vAffected[3]||vAffected[4]||vAffected[5]))
      return dGstar;

   // Direct-space metric tensor G and its derivative dG/dpar, using the
   // effective (constrained) lattice parameters
   const CrystVector_REAL cell=this->GetLatticePar();
   const REAL a=cell(0), b=cell(1), c=cell(2);
   const REAL ca=cos(cell(3)), cb=cos(cell(4)), cg=cos(cell(5));
   const REAL sa=sin(cell(3)), sb=sin(cell(4)), sg=sin(cell(5));
   const REAL G[3][3]={{a*a   , a*b*cg, a*c*cb},
                       {a*b*cg, b*b   , b*c*ca},
                       {a*c*cb, b*c*ca, c*c   }};
   REAL dG[3][3]={{0,0,0},{0,0,0},{0,0,0}};
   if(vAffected[0]) {dG[0][0]+=2*a; dG[0][1]+=b*cg; dG[1][0]+=b*cg; dG[0][2]+=c*cb; dG[2][0]+=c*cb;}
   if(vAffected[1]) {dG[1][1]+=2*b; dG[0][1]+=a*cg; dG[1][0]+=a*cg; dG[1][2]+=c*ca; dG[2][1]+=c*ca;}
   if(vAffected[2]) {dG[2][2]+=2*c; dG[0][2]+=a*cb; dG[2][0]+=a*cb; dG[1][2]+=b*ca; dG[2][1]+=b*ca;}
   if(vAffected[3]) {dG[1][2]+=-b*c*sa; dG[2][1]+=-b*c*sa;}
   if(vAffected[4]) {dG[0][2]+=-a*c*sb; dG[2][0]+=-a*c*sb;}
   if(vAffected[5]) {dG[0][1]+=-a*b*sg; dG[1][0]+=-a*b*sg;}

   // G* = G^-1, and dG*/dpar = -G^-1 (dG/dpar) G^-1
   REAL Ginv[3][3];
   {
      const REAL det= G[0][0]*(G[1][1]*G[2][2]-G[1][2]*G[2][1])
                     -G[0][1]*(G[1][0]*G[2][2]-G[1][2]*G[2][0])
                     +G[0][2]*(G[1][0]*G[2][1]-G[1][1]*G[2][0]);
      Ginv[0][0]= (G[1][1]*G[2][2]-G[1][2]*G[2][1])/det;
      Ginv[0][1]=-(G[0][1]*G[2][2]-G[0][2]*G[2][1])/det;
      Ginv[0][2]= (G[0][1]*G[1][2]-G[0][2]*G[1][1])/det;
      Ginv[1][0]=Ginv[0][1];
      Ginv[1][1]= (G[0][0]*G[2][2]-G[0][2]*G[2][0])/det;
      Ginv[1][2]=-(G[0][0]*G[1][2]-G[0][2]*G[1][0])/det;
      Ginv[2][0]=Ginv[0][2];
      Ginv[2][1]=Ginv[1][2];
      Ginv[2][2]= (G[0][0]*G[1][1]-G[0][1]*G[1][0])/det;
   }
   REAL tmp[3][3];// (dG/dpar) G^-1
   for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
         tmp[i][j]=dG[i][0]*Ginv[0][j]+dG[i][1]*Ginv[1][j]+dG[i][2]*Ginv[2][j];
   for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
         dGstar(i,j)=-(Ginv[i][0]*tmp[0][j]+Ginv[i][1]*tmp[1][j]+Ginv[i][2]*tmp[2][j]);
   return dGstar;
}

void UnitCell::Init(const REAL a, const REAL b, const REAL c, const REAL alpha,
                    const REAL beta, const REAL gamma,const string &SpaceGroupId,
                    const string& name)
{
   VFN_DEBUG_ENTRY("UnitCell::Init(a,b,c,alpha,beta,gamma,Sg,name)",10)
   //mSpaceGroup.Print();
   mSpaceGroup.ChangeSpaceGroup(SpaceGroupId);
   //mSpaceGroup.Print();
   mCellDim(0)=a;
   mCellDim(1)=b;
   mCellDim(2)=c;
   mCellDim(3)=alpha;
   mCellDim(4)=beta;
   mCellDim(5)=gamma;
   mClockLatticePar.Click();
   this->UpdateLatticePar();
   if((mCellDim(3)<=0)||(mCellDim(3)>=M_PI)) throw ObjCrystException("alpha must be within ]0;pi[");
   if((mCellDim(4)<=0)||(mCellDim(4)>=M_PI)) throw ObjCrystException("beta must be within ]0;pi[");
   if((mCellDim(5)<=0)||(mCellDim(5)>=M_PI)) throw ObjCrystException("gamma must be within ]0;pi[");

   mClockMetricMatrix.Reset();
   mClockLatticeParUpdate.Reset();

   this->InitRefParList();
   this->InitMatrices();
   this->SetName(name);

   VFN_DEBUG_EXIT("UnitCell::Init(a,b,c,alpha,beta,gamma,Sg,name):End",10)
}

void UnitCell::InitOptions()
{
   VFN_DEBUG_ENTRY("UnitCell::InitOptions",10)
   static string ConstrainLatticeToSpaceGroupName;
   static string ConstrainLatticeToSpaceGroupChoices[2];

   static bool needInitNames=true;
   if(true==needInitNames)
   {
      ConstrainLatticeToSpaceGroupName="Constrain Lattice to SpaceGroup Symmetry";
      ConstrainLatticeToSpaceGroupChoices[0]="Yes (Default)";
      ConstrainLatticeToSpaceGroupChoices[1]="No (Allow Crystallographic Pseudo-Symmetry)";

      needInitNames=false;//Only once for the class
   }
   VFN_DEBUG_MESSAGE("UnitCell::Init(a,b,c,alpha,beta,gamma,Sg,name):Init options",5)
   mConstrainLatticeToSpaceGroup.Init(2,&ConstrainLatticeToSpaceGroupName,
                                        ConstrainLatticeToSpaceGroupChoices);
   mConstrainLatticeToSpaceGroup.SetChoice(0);
   this->AddOption(&mConstrainLatticeToSpaceGroup);
   VFN_DEBUG_EXIT("UnitCell::InitOptions",10)
}

void UnitCell::InitMatrices() const
{
   //:NOTE: The Matrices must remain upper triangular, since this is assumed for
   //optimization purposes in some procedures.
   if(  (mClockMetricMatrix>mClockLatticePar)
      &&(mClockMetricMatrix>mSpaceGroup.GetClockSpaceGroup())) return;//no need to update
   //this->UpdateLatticePar(); we should be able to do this...

   VFN_DEBUG_MESSAGE("UnitCell::InitMatrices() for crystal : "+this->GetName(),5)
   //mClockMetricMatrix.Print();
   //mClockLatticePar.Print();

   REAL a,b,c,alpha,beta,gamma;//direct space parameters
   REAL aa,bb,cc,alphaa,betaa,gammaa;//reciprocal space parameters
   REAL v;//volume of the unit cell
   a=this->GetLatticePar(0);
   b=this->GetLatticePar(1);
   c=this->GetLatticePar(2);
   alpha=this->GetLatticePar(3);
   beta=this->GetLatticePar(4);
   gamma=this->GetLatticePar(5);

   //cout <<a<<" "<<b<<" "<<c<<" "<<alpha<<" "<<beta<<" "<<gamma<<endl;

   v=sqrt(1-cos(alpha)*cos(alpha)-cos(beta)*cos(beta)-cos(gamma)*cos(gamma)
            +2*cos(alpha)*cos(beta)*cos(gamma));

   aa=sin(alpha)/a/v;
   bb=sin(beta )/b/v;
   cc=sin(gamma)/c/v;

   alphaa=acos( (cos(beta )*cos(gamma)-cos(alpha))/sin(beta )/sin(gamma) );
   betaa =acos( (cos(alpha)*cos(gamma)-cos(beta ))/sin(alpha)/sin(gamma) );
   gammaa=acos( (cos(alpha)*cos(beta )-cos(gamma))/sin(alpha)/sin(beta ) );

   //cout <<aa<<" "<<bb<<" "<<cc<<" "<<alphaa<<" "<<betaa<<" "<<gammaa<<endl;

   mBMatrix = aa ,  bb*cos(gammaa) , cc*cos(betaa) ,
               0  , bb*sin(gammaa) ,-cc*sin(betaa)*cos(alpha),
               0  , 0              ,1/c;
   //cout <<"B Matrix :"<<endl<< mBMatrix <<endl;

   mOrthMatrix = a  , b*cos(gamma) , c*cos(beta) ,
                 0  , b*sin(gamma) ,-c*sin(beta)*cos(alphaa),
                 0  , 0              ,1/cc;

   mOrthMatrixInvert=InvertMatrix(mOrthMatrix);
   mBMatrixInvert=InvertMatrix(mBMatrix);
   //cout << "Orth Matrix :"<<endl<<mOrthMatrix <<endl;
   //cout << "InvOrth Matrix :"<<endl<<mOrthMatrixInvert <<endl;
   //cout << "Orth * InvOrth matrix :"<<endl<<product(mOrthMatrix,mOrthMatrixInvert) <<endl;
   //cout << "InvOrth * Orth Matrix :"<<endl<<product(mOrthMatrixInvert,mOrthMatrix) <<endl;
   mClockMetricMatrix.Click();
   VFN_DEBUG_MESSAGE("UnitCell::InitMatrices():End.",5)
}

void UnitCell::UpdateLatticePar()
{
   if(  (mClockLatticeParUpdate>mSpaceGroup.GetClockSpaceGroup())
      &&(mClockLatticeParUpdate>mClockLatticePar)) return;
   VFN_DEBUG_ENTRY("UnitCell::UpdateLatticePar().",3)

   int num = mSpaceGroup.GetSpaceGroupNumber();
   if((num <=2)||(mConstrainLatticeToSpaceGroup.GetChoice()!=0))
   {
      mClockLatticeParUpdate.Click();
      return;
   }
   if((num <=15) && (0==mSpaceGroup.GetUniqueAxis()))
   {
      mCellDim(4)=M_PI/2.;
      mCellDim(5)=M_PI/2.;
      mClockLatticeParUpdate.Click();
      return;
   }
   if((num <=15) && (1==mSpaceGroup.GetUniqueAxis()))
   {
      mCellDim(3)=M_PI/2.;
      mCellDim(5)=M_PI/2.;
      mClockLatticeParUpdate.Click();
      return;
   }
   if((num <=15) && (2==mSpaceGroup.GetUniqueAxis()))
   {
      mCellDim(3)=M_PI/2.;
      mCellDim(4)=M_PI/2.;
      mClockLatticeParUpdate.Click();
      return;
   }
   if(num <=74)
   {
      mCellDim(3)=M_PI/2.;
      mCellDim(4)=M_PI/2.;
      mCellDim(5)=M_PI/2.;
      mClockLatticeParUpdate.Click();
      return;
   }
   if(num <= 142)
   {
      mCellDim(3)=M_PI/2.;
      mCellDim(4)=M_PI/2.;
      mCellDim(5)=M_PI/2.;
      mCellDim(1) = mCellDim(0) ;
      mClockLatticeParUpdate.Click();
      return;
   }
   if(mSpaceGroup.GetExtension()=='R')
   {
      mCellDim(4) = mCellDim(3);
      mCellDim(5) = mCellDim(3);
      mCellDim(1) = mCellDim(0);
      mCellDim(2) = mCellDim(0);
      mClockLatticeParUpdate.Click();
      return;
   }
   if(num <= 194) //||(mSpaceGroup.GetExtension()=='H')
   {//Hexagonal axes, for hexagonal and non-rhomboedral trigonal cells
      mCellDim(3) = M_PI/2.;
      mCellDim(4) = M_PI/2.;
      mCellDim(5) = M_PI*2./3.;
      mCellDim(1) = mCellDim(0) ;
      mClockLatticeParUpdate.Click();
      return;
   }
   mCellDim(3)=M_PI/2.;
   mCellDim(4)=M_PI/2.;
   mCellDim(5)=M_PI/2.;
   mCellDim(1) = mCellDim(0) ;
   mCellDim(2) = mCellDim(0) ;
   mClockLatticeParUpdate.Click();
   VFN_DEBUG_EXIT("UnitCell::UpdateLatticePar().",3)
   return;
}

void UnitCell::InitRefParList()
{
   VFN_DEBUG_ENTRY("UnitCell::InitRefParList()",5)
   //this->ResetParList();
   int num = mSpaceGroup.GetSpaceGroupNumber();
   bool a=true;
   bool b=true;
   bool c=true;
   bool alpha=true;
   bool beta=true;
   bool gamma=true;
   if(mConstrainLatticeToSpaceGroup.GetChoice()==0)
   {
      if(num <=2)
      {
      }
      else if((num <=15) && (0==mSpaceGroup.GetUniqueAxis()))
      {
         beta=false;
         gamma=false;
      }
      else if((num <=15) && (1==mSpaceGroup.GetUniqueAxis()))
      {
         alpha=false;
         gamma=false;
      }
      else if((num <=15) && (2==mSpaceGroup.GetUniqueAxis()))
      {
         alpha=false;
         beta=false;
      }
      else if(num <=74)
      {
         alpha=false;
         beta=false;
         gamma=false;
      }
      else if(num <= 142)
      {
         b=false;
         alpha=false;
         beta=false;
         gamma=false;
      }
      else if(mSpaceGroup.GetExtension()=='R')
      {//Rhombohedral
         b=false;
         c=false;
         alpha=true;
         beta=false;
         gamma=false;
      }
      else if(num <= 194)
      {//Hexagonal axes, for hexagonal and non-rhomboedral trigonal cells
         b=false;
         alpha=false;
         beta=false;
         gamma=false;
      }
      else
      {
         b=false;
         c=false;
         alpha=false;
         beta=false;
         gamma=false;
      }
   }
   REAL *pLatPar=mCellDim.data();
   if(this->GetNbPar()==0)
   {//:KLUDGE:
      {
         RefinablePar tmp("a",pLatPar,1.,100.,
                           gpRefParTypeUnitCellLength,REFPAR_DERIV_STEP_ABSOLUTE,
                           true,true,a,false,1.0);
         tmp.SetDerivStep(1e-4);
         tmp.AssignClock(mClockLatticePar);
         this->AddPar(tmp);
      }
      {
         RefinablePar tmp("b",pLatPar+1,1.,100.,
                           gpRefParTypeUnitCellLength,REFPAR_DERIV_STEP_ABSOLUTE,
                           true,true,b,false,1.0);
         tmp.SetDerivStep(1e-4);
         tmp.AssignClock(mClockLatticePar);
         this->AddPar(tmp);
      }
      {
         RefinablePar tmp("c",pLatPar+2,1.,100.,
                           gpRefParTypeUnitCellLength,REFPAR_DERIV_STEP_ABSOLUTE,
                           true,true,c,false,1.0);
         tmp.SetDerivStep(1e-4);
         tmp.AssignClock(mClockLatticePar);
         this->AddPar(tmp);
      }
      {
         RefinablePar tmp("alpha",pLatPar+3,.5,3.,
                           gpRefParTypeUnitCellAngle,REFPAR_DERIV_STEP_ABSOLUTE,
                           true,true,alpha,false,RAD2DEG);
         tmp.SetDerivStep(1e-4);
         tmp.AssignClock(mClockLatticePar);
         this->AddPar(tmp);
      }
      {
         RefinablePar tmp("beta",pLatPar+4,.5,3.,
                           gpRefParTypeUnitCellAngle,REFPAR_DERIV_STEP_ABSOLUTE,
                           true,true,beta,false,RAD2DEG);
         tmp.SetDerivStep(1e-4);
         tmp.AssignClock(mClockLatticePar);
         this->AddPar(tmp);
      }
      {
         RefinablePar tmp("gamma",pLatPar+5,.5,3.,
                           gpRefParTypeUnitCellAngle,REFPAR_DERIV_STEP_ABSOLUTE,
                           true,true,gamma,false,RAD2DEG);
         tmp.SetDerivStep(1e-4);
         tmp.AssignClock(mClockLatticePar);
         this->AddPar(tmp);
      }
   }
   else
   {//Just Fix the 'used' status
      this->GetPar(pLatPar+0).SetIsUsed(a);
      this->GetPar(pLatPar+1).SetIsUsed(b);
      this->GetPar(pLatPar+2).SetIsUsed(c);
      this->GetPar(pLatPar+3).SetIsUsed(alpha);
      this->GetPar(pLatPar+4).SetIsUsed(beta);
      this->GetPar(pLatPar+5).SetIsUsed(gamma);
   }

   VFN_DEBUG_EXIT("UnitCell::InitRefParList():Finished",5)
}

}
