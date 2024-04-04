 
  
  //results for both TPC and iTPC
  Double_t padRadii[100];
  int innerPadRows=-1;
  int outerPadRows=-1; 
  int nPadRows=-1;

 const Double_t outerSectorRowPitch     =    2.0; // ;

  void InitPadRadii(){

    innerPadRows   =   13; // ;
    outerPadRows   =   32; // ;
    nPadRows=innerPadRows+outerPadRows;

    const Double_t firstOuterSectorPadRow  =  127.195; // ;
    Double_t RadInner[13] = {   60,  64.8,  69.6,  74.4,  79.2,  84,  88.8,  93.6,  98.8,   104, 109.2, 114.4, 119.6};
    for (Int_t i = 0; i < innerPadRows; i++) padRadii[i]=RadInner[i];
    
    for (Int_t i = 0; i < outerPadRows; i++)
       padRadii[i+innerPadRows]   = firstOuterSectorPadRow + outerSectorRowPitch*i;
      
     for (Int_t i = 0; i < nPadRows; i++)
      cout<<"pad #"<<i<<" r="<< padRadii[i]<<endl;
   }

//------------------------------------
   //coppied from itpPadPlanes.C
 
  void InitPadRadii_iTPC(){

    innerPadRows   =   40; // ;
    outerPadRows   =   32; // ;
    nPadRows=innerPadRows+outerPadRows;


     int NinnerRows = 40;
    // int` row.innerPadRows       = 40 ; //            
  //row.superInnerPadRows      =          3; //            
  //row.innerSectorPadPitch    =       0.50;
   float innerSectorPadLength   =       1.55; //  
  const Double_t firstOuterSectorPadRow  =  127.195; // ;
      
  //row.innerSectorPadWidth    = row.innerSectorPadPitch  - 0.05;
  double innerSectorRowPitch1   = innerSectorPadLength + 0.05;  
  //row.innerSectorEdge      =     51.905; // ;          
  //row.innerSectorPadPlaneZ   =    209.707; // ;                  
  Double_t RImax = 118.2 + innerSectorRowPitch1/2.; // 118.2 + 1.6/2 = 119 cm  JT
  Double_t RImin = RImax - innerPadRows*innerSectorRowPitch1; // 55.8 - 1.6/2 = 55 cm JT
  Double_t ROmin = 127.195;
  float firstPadRow      = RImin + innerSectorRowPitch1/2; 
  for (Int_t i = 0; i < innerPadRows; i++) {
    padRadii[i]     = firstPadRow + innerSectorRowPitch1*i;
  }
   
    
    for (Int_t i = 0; i < outerPadRows; i++)
       padRadii[i+innerPadRows]   = firstOuterSectorPadRow + outerSectorRowPitch*i;
      
     for (Int_t i = 0; i < nPadRows; i++)
      cout<<"pad #"<<i<<" r="<< padRadii[i]<<endl;
   }



   int PadFromR(double r){
     int i=nPadRows-1;
     while (r<padRadii[i]) i--;
     //cout<<i<<" "<<padRadii[i]<<endl;
     if (i== nPadRows-1) return i;
     if (i>=innerPadRows){
      if ((r-padRadii[i])<outerSectorRowPitch/2.) return i; else return i+1;
     } 
     //inner sector
     if ( (r-padRadii[i])<padRadii[i+1]-r) return i; else return i+1;   
   
   }


  
  /* TPC
  row.ioSectorSeparation	 =      7.595; // ;
  row.firstRowWidth	 =     29.145; // ;
  row.lastRowWidth	 =      95.81; // ;
  row.innerPadRows	 =         13; // ;
  row.innerPadRows48	 =          8; // ;
  row.innerPadRows52	 =          5; // ;
  row.superInnerPadRows	 =          3; // ;
  row.innerSectorPadWidth	 =      0.285; // ;
  row.innerSectorPadLength	 =       1.15; // ;
  row.innerSectorPadPitch	 =      0.335; // ;
  row.innerSectorRowPitch1	 =        4.8; // ;
  row.innerSectorRowPitch2	 =        5.2; // ;
  row.innerSectorEdge	 =     51.905; // ;
  row.innerSectorPadPlaneZ	 =    209.707; // ;
  for (Int_t i = 0; i < row.innerPadRows; i++) {
    row.innerPadsPerRow[i] = nPadsInner[i];
    row.innerRowRadii[i]   = RadInner[i];
  }
  row.firstPadRow	 =  row.innerRowRadii[0];   
  row.outerPadRows	 =         32; // ;
  row.padRows	         = row.innerPadRows + row.outerPadRows;
  row.outerSectorPadWidth	 =       0.62; // ;
  row.outerSectorPadLength	 =       1.95; // ;
  row.outerSectorPadPitch	 =       0.67; // ;
  row.outerSectorLength	 =         62; // ;
  row.superOuterPadRows	 =          1; // ;
  row.firstOuterSectorPadRow	 =    127.195; // ;
  row.outerSectorRowPitch	 =          2; // ;
  row.outerSectorEdge	 =    121.732; // ;
  row.outerSectorPadPlaneZ	 =    210.107; // ;
  for (Int_t i = 0; i < row.outerPadRows; i++) {
    row.outerPadsPerRow[i] = nPadsOuter[i];
    row.outerRowRadii[i]   = row.firstOuterSectorPadRow + row.outerSectorRowPitch*i;
  }
  row.lastOuterSectorPadRow	 =  row.outerRowRadii[row.outerPadRows-1];   
  tableSet->AddAt(&row);
  // ----------------- end of code ---------------
  return (TDataSet *)tableSet;
}*/

/* iTPC
#include "tables/St_itpcPadPlanes_Table.h"

TDataSet *CreateTable() { 
  // -----------------------------------------------------------------
  // db/.const/StarDb/Geometry/tpc/.itpcPadPlanes/itpcPadPlanes Allocated rows: 1  Used rows: 1  Row size: 1392 bytes
  //  Table: itpcPadPlanes_st[0]--> itpcPadPlanes_st[0]
  // ====================================================================
  // ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_itpcPadPlanes")) return 0;
  itpcPadPlanes_st row;
  St_itpcPadPlanes *tableSet = new St_itpcPadPlanes("itpcPadPlanes",1);
  /* 
     Jim Thomas 04/02/15 http://rnc.lbl.gov/~jhthomas/public/iTPC/PadPlane/PadRow2HVandWire.pdf
     pads = 2*(row + 24 - Int_t ((row - 1)/7)); // row = 1 - 40;

Pads are spaced 0.5 x 1.6 cm with a 0.5 mm gap.
Which means the actual copper size of the pad is 0.45 x 1.55 cm.

     Jim Thomas 05/31/16
     pads = 2*(row + 25 - Int_t ((row - 1)/7)); // row = 1 - 40;
  *
  Int_t NinnerRows = 40;
  Int_t nPadsInner[40] = { //J.Thomas, 05/31/2016
    52, 54, 56, 58, 60, 62, 62, 64, 66, 68,
    70, 72, 74, 74, 76, 78, 80, 82, 84, 86,
    86, 88, 90, 92, 94, 96, 98, 98,100,102,
   104,106,108,110,110,112,114,116,118,120};
//   for (Int_t i = 0; i < NinnerRows; i++) {
//     Int_t r = i + 1;
//     nPadsInner[i] = 2*(r + 24 - Int_t ((r - 1)/7));
//   }

  Int_t nPadsOuter[32]  = {    98 ,    100,    102,    104,    106,    106,    108,    110,    112,    112,
             114,    116,    118,    120,    122,    122,    124,    126,    128,    128,
             130,    132,    134,    136,    138,    138,    140,    142,    144,    144,    
             144,    144};
  memset(&row,0,tableSet->GetRowSize());
  row.innerPadRows       = NinnerRows; //            
  row.superInnerPadRows      =          3; //            
  row.innerSectorPadPitch    =       0.50;
  row.innerSectorPadLength   =       1.55; //         
  row.innerSectorPadWidth    = row.innerSectorPadPitch  - 0.05;
  row.innerSectorRowPitch1   = row.innerSectorPadLength + 0.05;  
  row.innerSectorEdge      =     51.905; // ;          
  row.innerSectorPadPlaneZ   =    209.707; // ;                  
  Double_t RImax = 118.2 + row.innerSectorRowPitch1/2; // 118.2 + 1.6/2 = 119 cm  JT
  Double_t RImin = RImax - NinnerRows*row.innerSectorRowPitch1; // 55.8 - 1.6/2 = 55 cm JT
  Double_t ROmin = 127.195;
  row.firstPadRow      = RImin + row.innerSectorRowPitch1/2;
  for (Int_t i = 0; i < row.innerPadRows; i++) {
    row.innerPadsPerRow[i]   = nPadsInner[i];
    row.innerRowRadii[i]     = row.firstPadRow + row.innerSectorRowPitch1*i;
  }
  row.firstRowWidth      =   row.innerPadsPerRow[0]*row.innerSectorPadPitch;
  row.outerPadRows       =      32; // ;
  row.padRows              = row.innerPadRows + row.outerPadRows;
  row.outerSectorPadWidth    =    0.62; // ;
  row.outerSectorPadLength   =    1.95; // ;
  row.outerSectorPadPitch    = row.outerSectorPadWidth + 0.05;     
  row.outerSectorLength      =      62; // ;
  row.superOuterPadRows      =       1; // ;
  row.firstOuterSectorPadRow = ROmin;
  row.outerSectorRowPitch    = row.outerSectorPadLength + 0.05;    // ;
  row.outerSectorEdge      = 121.732; // ;
  row.outerSectorPadPlaneZ   = 210.107; // ;
  for (Int_t i = 0; i < row.outerPadRows; i++) {
    row.outerPadsPerRow[i] = nPadsOuter[i];
    row.outerRowRadii[i]   = row.firstOuterSectorPadRow + row.outerSectorRowPitch*i;
  }
  row.lastRowWidth     = row.outerPadsPerRow[row.outerPadRows-1]*row.outerSectorPadPitch;
  row.ioSectorSeparation   = 
    (row.outerRowRadii[                 0] - row.outerSectorRowPitch/2) -
    (row.innerRowRadii[row.innerPadRows-1] + row.innerSectorRowPitch1/2);
  row.lastOuterSectorPadRow  =  row.outerRowRadii[row.outerPadRows-1];   
  tableSet->AddAt(&row);
  // ----------------- end of code ---------------
  return (TDataSet *)tableSet;
}
*/
