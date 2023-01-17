#include "Selector.h"
#include "Nano.h"
#include "TSystem.h"

TChain chain("mkselectChain");
 
void mkchain(const char *h1dir = 0)
{
   if (h1dir) {
      gSystem->Setenv("H1",h1dir);
   }
   chain.SetCacheSize(20*1024*1024);
   chain.Add("$H1/dstarmb.root");
   chain.Add("$H1/dstarp1a.root");
   chain.Add("$H1/dstarp1b.root");
   chain.Add("$H1/dstarp2.root");
}