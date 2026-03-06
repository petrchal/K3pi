#ifndef K3PILIB
#define K3PILIB


/*
#include "TNtuple.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
*/
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TStyle.h"

//#include "TROOT.h"
#include "TEfficiency.h"
#include <iostream>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>
#include <ROOT/RDF/RResultMap.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RVec.hxx>
#include <ROOT/RTrivialDS.hxx>

#include<optional>
#include <vector>
#include <iterator>
//using namespace ROOT;

  typedef struct{
    const TString fileName;
    const TString lable;
    const int isMc=0;
    const char * trigger=NULL;// separate trigger cut
  } TFileDescription;

//------------------------------------------------------------------------
//----   definitions of plots---------------------------------
//------------------------------------------------------------------------
  //exp - name of the variable or expression to be plotted from TTree
  //title of the resulting histogram
  //lo, hi - histogram limits 0,0 means automatic - takes long
  //NOTE!!!!: do not use 0,0 for range of continuous variables - the division gets wrong results!!!!
  typedef struct{
      const char* expr;
      const char* title;
      const char* axisTitle=NULL;
      float lo; float hi;
      //comma-separated list of cuts which are disabled before plotting
      //work only when 
      const char* cutMods=NULL; //cut modification when plotting this histogram

    } TPlotDef;
  //I could also add info how to deal with it ..normalize, or not.. 

//same for 2D plots
typedef struct{
      const char* x_expr; //variable on x-axis
      const char* y_expr;//variable on y-axis
      const char* title;
      const char* x_axisTitle=NULL;
      const char* y_axisTitle=NULL;
      float x_lo; float x_hi;
      float y_lo; float y_hi;
      //comma-separated list of cuts which are disabled before plotting
      //work only when 
      const char* cutMods=NULL; //cut modification when plotting this histogram

    } TPlotDef_2D;

  


  //list of all plot definitions 
  bool gIgnoreCutMods=false; //setting to true can speed running - reduces data filtering
  typedef std::vector<TPlotDef> TPlotDefinitions;
  typedef std::vector<TPlotDef_2D> TPlotDefinitions_2D;



//-----------------------------------------------------------------------
//---structures for storing results----------------------------------
//-----------------------------------------------------------------------

template <typename T> struct TSinglePlotRes{
     //ugly way how to go around non-existence of default constructor for RResultMap
     // only first member is filled
     // in c++17 one can use std::optional 
     ROOT::RDF::Experimental::RResultMap<T> &ResMap(){return resMap[0];}
     std::vector<ROOT::RDF::Experimental::RResultMap<T>> resMap;
     //list of labels which identifity plot
     //now it dataId(file), modifer, variation 
     //not fully used so far
     std::tuple<string,string,string> labels;
     //TString label; //how to label in the graph    
     unsigned int normalization=0;

};

//List of results coming from the same definition.
//The list is usualy over different files
template <typename T, typename defT> struct TResultStack{ 
          defT def; 
          std::vector<TSinglePlotRes<T>> singlePlots; 
          //std::vector<TString> localcut;
          THStack *stack=NULL; 
          THStack *ratios=NULL;
        }; 

//results from one definition over all files(data imputs)
typedef TResultStack<TH1D,TPlotDef> SingleDefResStack1D;
typedef TResultStack<TH2D,TPlotDef_2D> SingleDefResStack2D;

//list of results - basicaly vector holding information about current position of interator
template<typename T>
class TResVector : public std::vector<T> {
private:
    typename std::vector<T>::iterator currentPos; // Iterator to track position

public:
    // Default constructor
    TResVector() : std::vector<T>(), currentPos(this->begin()) {}
    // Constructor with initial size
    TResVector(size_t size) : std::vector<T>(size), currentPos(this->begin()) {}
    // Constructor with size and default value
    TResVector(size_t size, const T& value) : std::vector<T>(size, value), currentPos(this->begin()) {}
    // Copy constructor
    TResVector(const TResVector<T>& other) : std::vector<T>(other), currentPos(this->begin()) {}
    // Move constructor
    TResVector(TResVector<T>&& other) noexcept : std::vector<T>(std::move(other)), currentPos(this->begin()) {}
    // Constructor from std::vector
    TResVector(const std::vector<T>& other) : std::vector<T>(other), currentPos(this->begin()) {}
    // Move constructor from std::vector
    TResVector(std::vector<T>&& other) noexcept : std::vector<T>(std::move(other)), currentPos(this->begin()) {}
    
    // Method to reset the position
    void resetPosition() {currentPos = this->begin();
    }

    // Method to advance the position by n steps
    void advancePosition(size_t n) {
        if (currentPos + n <= this->end()) {
            currentPos += n;
        } else {
            currentPos = this->end();
        }
    }

    // Method to get the current position
    typename std::vector<T>::iterator *CurrentPosition() {
        return &currentPos;
    }

    // Method to print the current position
    void printCurrentPosition() const {
        if (currentPos != this->end()) {
            std::cout << "Current position points to position: - not end"  << std::endl;
        } else {
            std::cout << "Current position is at the end of the vector." << std::endl;
        }
    }
  
};

//list over all plot definitions
typedef TResVector<SingleDefResStack1D> ResultList1D; 
typedef TResVector<SingleDefResStack2D> ResultList2D; 


//----------------------------------------------------
// for results on efficiciency
//----------------------------------------------------


 // structure contains results of one efficiencey plot for single variation from all files
   //.. Note:the THStack is useful only for 1D histos 
   template <typename T> struct TEffResult{

     //vector over files (data sets) of plots     
       std::vector<TEfficiency *>eff;
       std::vector<TH1D *>      hNum; 
       std::vector<TH1D *>      hDen; 
       std::vector<TH1D *>      ratio; //efficiency other way

   };

    
/*
   // efficiencies from one definition from multiple files (data sets)
   template <typename T> struct TEffFromSingleDef{ 
          TPlotDef def; 

         //workaround - num and den belong to TEFFResults
       
          std::vector<TEffResult<T>> res;

          //ugly way how to go around non-existence of default constructor for RResultMap
          //vector over files
          std::vector<ROOT::RDF::Experimental::RResultMap<T>> num;
          std::vector<ROOT::RDF::Experimental::RResultMap<T>> den;
       
          /*
          THStack *stack_num=NULL; 
          THStack *stack_den=NULL;
          THStack *stack_ratio=NULL;
         * /

          THStack stack_num; 
          THStack stack_den;
          THStack stack_ratio;
        }; 


typedef std::vector<TEffFromSingleDef<TH1D>> TEffList;
*/


// efficiencies from one definition vectored over all variations (and files)
template <typename T> struct TEffFromSingleDef{ 
          //------- part filled by AddPlot ----------------
          TPlotDef def; 

          //ugly way how to go around non-existence of default constructor for RResultMap
          //vector over files
          std::vector<ROOT::RDF::Experimental::RResultMap<T>> numMap;
          std::vector<ROOT::RDF::Experimental::RResultMap<T>> denMap;
         // -------------------------------------------------------

         // --- part of results filled  during plotting 

          //vectors over variations      
          std::vector<TEffResult<T>*> res;

        
          /*
          THStack *stack_num=NULL; 
          THStack *stack_den=NULL;
          THStack *stack_ratio=NULL;
         */

          //results combined for plotting
          THStack stack_num; 
          THStack stack_den;
          THStack stack_ratio;
        }; 
  
typedef std::vector<TEffFromSingleDef<TH1D>> TEffList;

//-----------------------------------------------------------------------
//---structures for storing CUTS and cut variations----------------------
//-----------------------------------------------------------------------

class K3PiCut{
public:
   
   void Clear(){m.clear();}
   string Str();//returns string of the compiled cut;
   string operator[]( char *v) {return m[v];}
   string& operator [] (const char* i) {return m[string(i)];}
  /*int H[100];
    int operator [] (int i) const {return H[i];}
    int& operator [] (int i) {return H[i];}
    */
    K3PiCut operator + (const K3PiCut& obj);
    K3PiCut operator - (const K3PiCut& obj);
    int Remove(const char* var);
    int Replace(const char* var);
 // private:  
   std::map<std::string,std::string> m;
  };

 K3PiCut K3PiCut::operator + (const K3PiCut& obj){
  K3PiCut res=*this;
   res.m.insert(obj.m.begin(),obj.m.end());
  return res;
 }
   
 K3PiCut K3PiCut::operator - (const K3PiCut& obj){
  K3PiCut res=*this;
  for (auto o:obj.m) res.m.erase(o.first);
  return res;
 }
  
//remove explicitly cut by name - returns how many found
//takes comma separated list;
 int K3PiCut::Remove(const char* var){
    if (var==NULL) return 0;
    std::stringstream test(var);
    std::string segment; std::vector<std::string> seglist;
    while(std::getline(test, segment, ',')) seglist.push_back(segment);
    int n=0;
    for (auto s:seglist) n+=m.erase(s.c_str());
    return n;
  }
  
//replace cut by name - returns how many found
//takes comma separated list name:cut
int K3PiCut::Replace(const char* var){
    if (var==NULL) return 0;
    std::stringstream test(var);
    std::string segment; std::vector<std::string> seglist;
    while(std::getline(test, segment, ',')) seglist.push_back(segment);
    int n=0;
    for (auto s:seglist){
        int pos=s.find(":");
        if (pos<0) continue;
        std::string name = s.substr(0,pos);
        s.erase(0, pos + 1);
        n+=m.erase(name.c_str()); 
        m[name.c_str()]=s.c_str();
    } 
    return n;
  }

  //put all cuts together
  std::string K3PiCut::Str(){
         string res="";
         if (m.size()==0) return res;
         for( auto it = m.begin(); it != m.end(); ++it )
           {
           if (it!=m.begin()) res+="&&";
           //cout << it->first; // key
           res += it->second;
           }
         return res;
  }

  
//------cut variations-------------------------------------------------------
//https://root.cern/doc/master/classROOT_1_1RDF_1_1RInterface.html#a9b67e8eb7

//---cut variables for varying
class TNewVariables: public K3PiCut{
public:
  ROOT::RDF::RNode DefineNewVariables(ROOT::RDF::RNode nod);
};

ROOT::RDF::RNode TNewVariables::DefineNewVariables(ROOT::RDF::RNode nod){
  for (auto v : m) {
    std::cout << "Defining new variable: " << v.first << "=" << v.second << std::endl;
    try {
        nod = nod.Define(v.first, v.second);
    } catch (const std::exception &e) {
        std::cerr << "Error defining variable " << v.first
                  << ": " << e.what() << std::endl;
    }
}
return nod;
}



TNewVariables NewVars;


ROOT::RDF::RNode DefineNewVariables(ROOT::RDF::RNode nod){
   return NewVars.DefineNewVariables(nod);
}

//---------------------------
typedef struct{
    const char*  var;
    const std::vector<string> variations;
    } TCutVariation;

const int cMaxVariationStyles=4;
int cVariationStyle[cMaxVariationStyles]={kSolid,kDashed,kDotted,kDashDotted};
int GetVariationLineStyle(int i){ if (i>=cMaxVariationStyles) return kDashed; return cVariationStyle[i];}

//--------------------------
ROOT::RDF::RNode AddVariations(TCutVariation& varyDefs, ROOT::RDF::RNode node)
{ //df.Vary({"x", "y"}, "ROOT::RVec<ROOT::RVecD>{{x*0.9, x*1.1}, {y*0.9, y*1.1}}", 2, "xy")
   cout<<"AddVariations:"<<endl;
  //  TString vec="ROOT::RVec<ROOT::RVecD>{{"; from manual
  TString vec="ROOT::RVecD{";
 for (const auto& value : varyDefs.variations){vec+=value;vec+=",";}
  int len = vec.Length(); vec.Remove(len-1,1); //remove last ",""
  vec+="}";
  cout<<"  "<<vec<<endl;
  auto  res=node.Vary(varyDefs.var,vec.Data(),varyDefs.variations,varyDefs.var);
    //auto  res=node.Vary(varyDefs.var,vec.Data(),varyDefs.variations.size(),varyDefs.var);
  return res;
}

//======================================================
//===========plotting==================================
//======================================================


//add plot defintion to the list that is to be processed
void AddPlots_1D(TPlotDefinitions& plotDefs, ROOT::RDF::RNode node,ResultList1D& rlist,
 ResultList1D::iterator &iter, const char * prefix = "", bool rebin=false,bool ignoreRange=false)
 {
    int tmpVarCount=0;
      for (auto def: plotDefs){ //loop over plot definiton and book the plots
        //if collumn name existed I woudl not need to do the Define.. but how to simply find out?
        TString var="tmpVar";var+=tmpVarCount++;
        TString title=prefix; title+=": "; title+=def.title;
        auto h=node.Define(var.Data(),def.expr) //this must be done for calculated variables
                   .Histo1D(ROOT::RDF::TH1DModel(def.expr,title, 100/rebin, (ignoreRange)?0:def.lo, (ignoreRange)?0:def.hi),var.Data()); 

        //If at the end of list add Results stack for another variable          
        if (iter==rlist.end()){ 
            SingleDefResStack1D r;
            r.def=def; r.stack=NULL;
            cout<<"pushing "<<def.expr<<endl;
            rlist.push_back(r);
            iter=rlist.end(); --iter;
            //std::advance(iter, rlist.size());
        }
        //iter is now pointing to valid TResultStack
        //add result from current file(node) for current variable
        TSinglePlotRes<TH1D> sr;
        sr.resMap.push_back(ROOT::RDF::Experimental::VariationsFor(h)); 
        sr.labels={prefix,"ss",""};
        //cout<<" currently at "<<it1D->def.expr<<" "<<it1D-rlist.begin()<<endl;
        iter++; //increase position int list of plotted variables
  }
}

//-----------------------------------------------------------------------
//this not only adds plots, but for each plots modifies the cut so that it can be plotted without bounds
ROOT::RDF::RNode AddPlots4QA(TPlotDefinitions& plotDefs, ROOT::RDF::RNode node,K3PiCut defaultCut, ResultList1D& rlist,
 //ResultList1D::iterator &iter, 
 const char * plotprefix = "", const char * fileprefix = "",bool rebin=false,bool ignoreRange=false, unsigned long normalization=0)
 {
   
   auto iter=rlist.CurrentPosition();
   //rlist->printCurrentPosition();

    ROOT::RDF::RNode defaultNode=node.Filter(defaultCut.Str()); //nod with all cuts applied
    int tmpVarCount=0;
   

    for (auto def: plotDefs){ //loop over single plot definiton and book the plots
         //for first file this means to create the structure, others only add histogram
     
         if (*iter==rlist.end()){ 
            SingleDefResStack1D r;
            r.def=def; r.stack=NULL;
            cout<<"pushing "<<def.expr<<endl;
            rlist.push_back(r);
            *iter=rlist.end();(*iter)--; //last element
          }
         //iter is now poiting to valid TResultStack
        
        
        //if collumn name existed I woudl not need to do the Define.. but how to simply find out?
        //column names not working in recent ROOT version - maybe later
        TString var="tmpVar";var+=tmpVarCount++;
        TString title=def.title; title+=" - "; title+=plotprefix; 
        auto h=defaultNode.Define(var.Data(),def.expr) //this must be done for calculated variables 
                   .Histo1D(ROOT::RDF::TH1DModel(def.expr,title, 100/rebin, (ignoreRange)?0:def.lo, (ignoreRange)?0:def.hi),var.Data()); 
       //now add it to the results at current position
        TSinglePlotRes<TH1D> sr;
        sr.resMap.push_back(ROOT::RDF::Experimental::VariationsFor(h)); 
        sr.labels={fileprefix,plotprefix,""};
        sr.normalization=normalization;
       
        //add second one for comparison with a given cut
        
        if (def.cutMods && (gIgnoreCutMods==false)){
          TSinglePlotRes<TH1D> sr;
          K3PiCut cut=defaultCut;
          cout<<"  modified:"<<def.cutMods<<endl;
          cout<<"        "<<cut.Str()<<endl;
          cout<<cut.Remove(def.cutMods)<<endl;
          cout<<"        "<<cut.Str()<<endl<<endl;
          auto h2=node.Filter(cut.Str()).Define(var.Data(),def.expr) //with modified cut
                  .Histo1D(ROOT::RDF::TH1DModel(def.expr,title, 100/rebin, (ignoreRange)?0:def.lo, (ignoreRange)?0:def.hi),var.Data()); 
          sr.resMap.push_back(ROOT::RDF::Experimental::VariationsFor(h2));
          TString bb;//=fileprefix;
          bb+="disabled ";bb+=def.cutMods; bb+=" cut";
         // sr.label=bb; 
          sr.labels={fileprefix,plotprefix,bb.Data()};
          sr.normalization=normalization;
      
          (*iter)->singlePlots.push_back(sr); //use prefix as a label of the histogram
        } 
        
        //because of drawing order push first the version with cut modification
        (*iter)->singlePlots.push_back(sr); //use prefix as a label of the histogram
       
     //  if (def.cutMods && iter->singlePlots.size()>1) goto SKIP;
     //the default one with full cut, but only for fist file
   
     //SKIP:  
       (*iter)++; //increase position int list of plotted variables
      }

      return defaultNode;
}

//-----version for 2D plots--------------------------------
//this not only adds plots, but for each plots modifies the cut so that it can be plotted without bounds
ROOT::RDF::RNode AddPlots4QA(TPlotDefinitions_2D& plotDefs, ROOT::RDF::RNode node,K3PiCut defaultCut, ResultList2D& rlist,
 //ResultList2D::iterator &iter, 
 const char * plotprefix = "", const char * fileprefix = "",bool rebin=false,bool ignoreRange=false,unsigned int normalization=0)
 {
    auto iter=rlist.CurrentPosition();

    ROOT::RDF::RNode defaultNode=node.Filter(defaultCut.Str()); //nod with all cuts applied
    int tmpVarCount=0;
    
    for (auto def: plotDefs){ //loop over single plot definiton and book the plots
         //for first file this means to create the structure, others only add histogram
     
         if (*iter==rlist.end()){ 
            SingleDefResStack2D r;
            r.def=def; r.stack=NULL;
            cout<<"pushing 2D plot of "<<def.x_expr<<" : "<<def.y_expr<<endl;
            rlist.push_back(r);
            (*iter)=rlist.end();(*iter)--; //last element
          }
         //iter is now poiting to valid TResultStack
        
        
        TString var="tmpVar";var+=tmpVarCount++;
        TString var2="tmpVar";var2+=tmpVarCount++;
        TString exp2D=def.x_expr;exp2D+=":";exp2D+=def.y_expr;
        TString title=def.title; title+=" - "; title+=plotprefix; 
        //this gets realy complicated for 2D 
        auto h=defaultNode.Define(var.Data(),def.x_expr).Define(var2.Data(),def.y_expr)  //this must be done for calculated variables 
                   .Histo2D(ROOT::RDF::TH2DModel(exp2D.Data(),title, 
                     100/rebin, (ignoreRange)?0:def.x_lo, (ignoreRange)?0:def.x_hi, 100/rebin, (ignoreRange)?0:def.y_lo, (ignoreRange)?0:def.y_hi),
                     var.Data(),var2.Data()); 
       //now add it to the results at current position
        TSinglePlotRes<TH2D> sr;
        sr.resMap.push_back(ROOT::RDF::Experimental::VariationsFor(h)); 
        sr.labels={fileprefix,plotprefix,""};
        sr.normalization=normalization;
        //add second one for comparison with a given cut
        
        if (def.cutMods && (gIgnoreCutMods==false)){
          TSinglePlotRes<TH2D> sr;
          K3PiCut cut=defaultCut;
          cout<<"  modified:"<<def.cutMods<<endl;
          cout<<"        "<<cut.Str()<<endl;
          cout<<cut.Remove(def.cutMods)<<endl;
          cout<<"        "<<cut.Str()<<endl<<endl;
          //with modified cut
          auto h2=node.Filter(cut.Str()). 
                Define(var.Data(),def.x_expr).Define(var2.Data(),def.y_expr)  //this must be done for calculated variables 
                   .Histo2D(ROOT::RDF::TH2DModel(exp2D.Data(),title, 
                     100/rebin, (ignoreRange)?0:def.x_lo, (ignoreRange)?0:def.x_hi, 100/rebin, (ignoreRange)?0:def.y_lo, (ignoreRange)?0:def.y_hi),
                     var.Data(),var2.Data()); 

          sr.resMap.push_back(ROOT::RDF::Experimental::VariationsFor(h2));
          TString bb;//=fileprefix;
          bb+="disabled ";bb+=def.cutMods; bb+=" cut";
         // sr.label=bb; 
          sr.labels={fileprefix,plotprefix,bb.Data()};
          sr.normalization=normalization;
          (*iter)->singlePlots.push_back(sr); //use prefix as a label of the histogram
        } 
        
        //because of drawing order push first the version with cut modification
        (*iter)->singlePlots.push_back(sr); //use prefix as a label of the histogram
       
     //  if (def.cutMods && iter->singlePlots.size()>1) goto SKIP;
     //the default one with full cut, but only for fist file
   
     //SKIP:  
       (*iter)++; // position int list of plotted variables
      }

      return defaultNode;
}



//-----------------------------------------------------------------------
//histograms and normalized ratios for one variable defintion over files ( and variations)
void Draw1Dstack(SingleDefResStack1D &r,const char* variation="nominal", const bool renormalize=false){

  const int   color[]={kBlack,kBlue,kRed,kGreen,kMagenta,kCyan}; 


    //auto hist=r.singlePlots.begin()->hist->GetPtr();
    //r.stack=new  THStack(hist); //copy settign from first histogram ...notw
    //auto s=new  THStack(const TH1* hist, Option_t* axis = "x", const char* name = 0, const char* title = 0, Int_t firstbin = 1, Int_t lastbin = -1, Int_t firstbin2 = 1, Int_t lastbin2 = -1, Option_t* proj_option = "", Option_t* draw_option = "");
  
    r.stack=new THStack();//r.def.expr,r.def.title);
    r.stack->SetTitle(r.def.expr);
      
    r.ratios=new THStack("","");

    auto c=new TCanvas("","",800,600);


    Double_t pdiv = 0.3;
    TPad *pad1 = new TPad("p1", "p1", 0., pdiv, 1., 1.); // upper
    TPad *pad2 = new TPad("p2", "p2", 0., 0., 1., pdiv); // lower
    pad1->Draw();
    pad2->Draw();

    pad1->cd();

  //global constants with range of the plot
  //const Double_t gxmin = 17., gymin = 15., gxmax = 1990., gymax = 1000.;

  // TH1F* frame1 = gPad->DrawFrame(gxmin,gymin,gxmax,gymax);

 //  frame1->GetXaxis()->SetMoreLogLabels();
    gPad->SetLeftMargin(0.09);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.11);//0.025
    gPad->SetBottomMargin(0.);

    pad2->SetLeftMargin(0.09);
    pad2->SetRightMargin(0.05);
    pad2->SetBottomMargin(0.2);//0.025
    pad2->SetTopMargin(0.);


    pad1->cd();
    auto l = new TLegend(0.65,0.75,0.9,0.9);
    l->SetHeader("","C");

    int id=0;//ugly but...


    TH1* h_den;
    bool first=true;
    for ( auto &sr : r.singlePlots){ //over files- TSinglePlotRes
      auto  keys=sr.ResMap().GetKeys();
      int  nVariations =keys.size();

      TString str=get<0>(sr.labels);str+=" ";str+=get<2>(sr.labels);
      //cout<<"sr.labels="<<str.Data()<<endl;
      l->AddEntry((TObject*)0, "", "");  
      if (nVariations>1) l->AddEntry((TObject*)0, str.Data(), "");
     
      for (int iv=0;iv<nVariations;iv++){ //loop over variations
         cout<<"  variation "<<keys[iv]<<endl;
         pad1->cd();
         if (renormalize && (sr.normalization>0)){ //forced to renormalize the plots - usually by number of events
            cout<<" !!!! renormalizing histograms per user request!!!"<<endl;
            sr.ResMap()[keys[iv]].Scale(1./((Double_t)sr.normalization));
         }
         TH1* h_copy=(TH1*)sr.ResMap()[keys[iv]].Clone(); //use default variation hist[]

         //if (iv==0){ //this will be denumerator
         if (first){ 
           cout<<" making denum"<<endl;
           h_den=(TH1*)h_copy->Clone();
           h_den->Scale(1./h_den->Integral());
         } 

        //rebin and scale
        // if(rebin>1) h_copy->Rebin(rebin);
        // if (normalize) h_copy->Scale(1./h_copy->Integral());
        //color and other things
        h_copy->SetTitle(r.def.expr);
        //h_copy->SetLineColor(color[id]);
        h_copy->SetLineColor(id+1);
        h_copy->SetMarkerColor(id+1);
       //add to stack and label
        r.stack->Add(h_copy);
        //l->AddEntry(h_copy,get<0>(sr.labels).c_str(),"l");
        if (nVariations>1) {
          TString lb=/*glb+=":: ";lb+=*/keys[iv].c_str();
          l->AddEntry(h_copy,lb,"l");
          }
          else  l->AddEntry(h_copy,str.Data(),"l");

      //normalized ratio
      //TSinglePlotRes<TH1D>* first=&(r.singlePlots.front());
      //if (&sr!=first){
     // if (iv>0){
      if (!first){
          cout<<" taking ratio"<<endl;
          TH1* h_num=(TH1*)h_copy->Clone();
          h_num->Scale(1./h_num->Integral());
          h_num->Divide(h_den);
          pad2->cd();
          r.ratios->Add(h_num);
        //h_num->Draw();
          } else cout<<" no ratio"<<endl;
         first=false; 
        id++;
      }

  } //iv loop 
  pad1->cd();

  r.stack->Draw("ehistnostack");
  r.stack->SetMinimum(0.001);

   Float_t siz = 0.045;
  
  auto frame1 = r.stack->GetHistogram();
  if (frame1) {
    frame1->SetTitleSize(siz);       frame1->SetLabelSize(siz);
    frame1->SetTitleSize(siz, "Y");  frame1->SetLabelSize(siz, "Y");
    frame1->SetLabelSize(0.00001, "X");
    frame1->GetXaxis()->SetTitleOffset(1.4);
    frame1->SetTitle(r.def.expr);
    
    
    l->Draw("same");
    //pad1->Modified();

    //bottom pad
    pad2->cd();
    r.ratios->Draw("enostack");
    if (r.ratios->GetXaxis()) r.ratios->GetXaxis()->SetTitle(r.def.axisTitle);
   }

  auto frame2 = r.ratios->GetHistogram();
  if (frame2) {
    //cout<<"min="<<r.ratios->GetMinimum()<<endl;
    ///cout<<"max="<<r.ratios->GetMaximum()<<endl;
    // cout<<"min2="<<r.ratios->GetYaxis()->GetXmin()<<endl;
    //cout<<"max2="<<r.ratios->GetYaxis()->GetXmax()<<endl;
    //if (r.ratios->GetYaxis()->GetXmin()<0) 

    //r.ratios->SetMinimum(0.5);
    //r.ratios->SetMaximum(2);

    siz = siz*(1.+(1.-pdiv));
    frame2->SetTitleSize(siz);       frame2->SetLabelSize(siz);
    frame2->SetTitleSize(siz, "Y");  frame2->SetLabelSize(siz, "Y");
    frame2->GetXaxis()->SetTitleOffset(1.1);
    frame2->GetYaxis()->SetTitle("ratio of normalized  ");
    frame2->GetYaxis()->SetTitleOffset(0.54); //0.64
   } 
    //pad2->Modified();

    c->Modified();
    c->Write(r.def.title);
  }

//---- simplifieed drawing for 2D histogram - so far does not handle comparision
void Draw2Dstack(SingleDefResStack2D &r,const char* variation="nominal"){

  const int   color[]={kBlack,kBlue,kRed,kGreen,kMagenta,kCyan}; 


    //auto hist=r.singlePlots.begin()->hist->GetPtr();
    //r.stack=new  THStack(hist); //copy settign from first histogram ...notw
    //auto s=new  THStack(const TH1* hist, Option_t* axis = "x", const char* name = 0, const char* title = 0, Int_t firstbin = 1, Int_t lastbin = -1, Int_t firstbin2 = 1, Int_t lastbin2 = -1, Option_t* proj_option = "", Option_t* draw_option = "");
      
    auto tmp = r.singlePlots.front().ResMap().GetKeys().size();
    int nn=r.singlePlots.size()*tmp;
    int nx=sqrt(nn); int ny=nx; 
    if (nx*ny < nn) nx++;
    if (nx*ny < nn) ny++;
    

    auto c=new TCanvas("","",nx*800,ny*600);
    c->Divide(nx,ny);

/*
 //  frame1->GetXaxis()->SetMoreLogLabels();
    gPad->SetLeftMargin(0.09);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.11);//0.025
    gPad->SetBottomMargin(0.);

    pad2->SetLeftMargin(0.09);
    pad2->SetRightMargin(0.05);
    pad2->SetBottomMargin(0.2);//0.025
    pad2->SetTopMargin(0.);
  */

    /*
    pad1->cd();
    auto l = new TLegend(0.65,0.75,0.9,0.9);
    l->SetHeader("","C");
        int id=0;//ugly but...
    */

      int pad=1;

     for ( auto &sr : r.singlePlots){ //over files- TSinglePlotRes
      TString str=get<0>(sr.labels);str+=" ";str+=get<2>(sr.labels);
      cout<<"sr.labels="<<str.Data()<<endl;
      //l->AddEntry((TObject*)0, "", "");  
     // l->AddEntry((TObject*)0, str.Data(), "");
      //denominator - te first of variation
    
      //loop over variations
      auto  keys=sr.ResMap().GetKeys();
      int  nVariations =keys.size();
      for (int iv=0;iv<nVariations;iv++){ 
      cout<<"  variation "<<keys[iv]<<endl;
       c->cd(pad++);


       TH1* h_copy=(TH1*)sr.ResMap()[keys[iv]].Clone(); //use default variation hist[]

        //color and other things
        TString titl=str.Data();
        if (iv>0) {titl+=" : ";titl+= keys[iv];}//r.def.x_expr;titl+=" vs ";titl+=r.def.y_expr;
        h_copy->SetTitle(titl);
        h_copy->GetXaxis()->SetTitle(r.def.x_axisTitle);
        h_copy->GetYaxis()->SetTitle(r.def.y_axisTitle);
         //add to stack and label
        //l->AddEntry(h_copy,get<0>(sr.labels).c_str(),"l");
        //TString lb=keys[iv].c_str();
        //l->AddEntry(h_copy,lb,"l");
        //h_copy->SetTitle(get<0>(sr.labels).c_str());
        h_copy->Draw("colz");
        gPad->SetLogz();
        gPad->Modified();
      } //over varitions

  } //over files loop

    c->Modified();
    c->Write(r.def.title);
  }

//-----------------------------------------------------------------------
 void DrawResults(ResultList1D &results, const bool renormalize=false, const char* whichVariation=NULL){
    cout<<endl<<"DrawResults"<<endl;
//draw 1D histograms in stack
    gStyle->SetPadTickY(1);
    gStyle->SetTickLength(0.02,"Y");
//gStyle->SetOptStat(0);

   for ( auto &r : results){ //over plots TResultStack, r-is for single plot definition
        cout<<"Plotting Stack for expression: "<<r.def.expr<<endl;
        if (r.singlePlots.size()==0){ cout<<"  !!!!plot empty"<<endl; continue; }
        //if (r.singlePlots.size()>1) 
          Draw1Dstack(r,whichVariation,renormalize);\
          /*
  else {
    cout<<"simple plot"<<endl;
    auto c=new TCanvas("","",800,600);
    gPad->Modified(); gPad->Update(); 
  // gSystem->ProcessEvents();
    auto keys=r.singlePlots[0].resMap[0].GetKeys();
    int  nVariations=keys.size();
    cout<<"variations="<<nVariations<<endl;
    for (int v=0;v<nVariations;v++){
      cout<<"   "<<keys[v]<<endl;
      TH1* h=(TH1*) r.singlePlots[0].ResMap()[keys[v]].Clone();
      h->SetDirectory(0);
      TString title=r.singlePlots.begin()->label; title+=": "; title+=r.def.title;
      h->SetTitle(title);
      h->GetXaxis()->SetTitle(r.def.expr);
    //c->cd();
    //h->DrawClone();
      if (v==0) h->Draw();else h->Draw("same"); 
    }
  }
  */
    }
}

 

//-----2D version--------------------------------
//sofar I cannot do the vgood plotting for variation of 2D histograms
  void DrawResults(ResultList2D &results,const char* whichVariation=NULL){
    cout<<endl<<"DrawResults"<<endl;
//draw 1D histograms in stack
    gStyle->SetPadTickY(1);
    gStyle->SetTickLength(0.02,"Y");
//gStyle->SetOptStat(0);

   for ( auto &r : results){ //over plots TResultStack, r-is for single plot definition
        cout<<"Plotting 2D Stack for "<<r.def.x_expr<<" : "<< r.def.y_expr<<endl;
        if (r.singlePlots.size()==0){ cout<<"  !!!!plot empty"<<endl; continue; }
        //if (r.singlePlots.size()>1) 
          Draw2Dstack(r,whichVariation);\
    }
}



//---------------------------------------------------
//--- code for efficiency plotting ------------------
//---------------------------------------------------




/*
//---progress bar ---------------------------------------------
 const UInt_t barWidth = 6;//60;
 ULong64_t processed = 0, totalEvents = 0;
 std::string progressBar;
 std::mutex barMutex;
 auto registerEvents = [](ULong64_t nIncrement) {totalEvents += nIncrement;};

 ROOT::RDF::RResultPtr<ULong64_t> AddProgressBar(ROOT::RDF::RNode df, int everyN=1000000, int totalN=10000000) {
  registerEvents(totalN);
  auto c = df.Count();
  c.OnPartialResultSlot(everyN, [everyN] (unsigned int slot, ULong64_t &cnt){
    std::lock_guard<std::mutex> l(barMutex);
            processed += everyN; //everyN captured by value for this lambda
            progressBar = "[";
            for(UInt_t i = 0; i < static_cast<UInt_t>(static_cast<Float_t>(processed)/totalEvents*barWidth); ++i){
              progressBar.push_back('|');
            }
            // escape the '\' when defined in python string
            std::cout << "\r" << std::left << std::setw(barWidth) << progressBar << "] " << processed << "/" << totalEvents << std::flush;
          });
  return c;
}
*/


#endif