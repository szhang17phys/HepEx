/// \file
/// \ingroup tutorial_fit
/// \notebook
/// Replacement test
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \author Yanwen Liu
#include <algorithm>
#include <TH1.h>

using namespace std;

int yw_marktwain() {

  // Input data
  const int ndata = 18;
  const int n1=8;
  const int n2=10;

  const double fracs[18] = {
    0.225, 0.262, 0.217, 0.240, 0.230, 0.229, 0.235, 0.217,
    0.209,0.205,0.196,0.210, 0.202,0.207,0.224, 0.223, 0.220, 0.201
  };

  // For replacement
  std::vector<bool> flag(ndata);
  std::fill(flag.begin(), flag.end(), false);
  std::fill(flag.end()-n2, flag.end(), true);

  // Fill a histogram for the test statistics
  TH1F *hDiff = new TH1F("hDiff", "diffs", 1000, -0.1, 0.1);
  int irun=0;
  double diff0 = 0.;
  int ntrials = 0;
  int np=0;

  // permutation: from 18 elements, sample 8
  do{
    double sum1=0, sum2=0;
    // Based on the flag to decide if the element is picked up or not: for a given element, it is either picked up or not
    for(int i=0; i<ndata; i++) {
      if(flag[i]) sum2 += fracs[i];
      else        sum1 += fracs[i];
    }
    double diff = sum1/float(n1) - sum2/float(n2);
    if(irun==0) {
      diff0 = diff;
      cout << diff0 << endl;
      cout << sum1/n1 << " " << sum2/n2 << endl;
    }
    // Check the test stat
    if(fabs(diff) > fabs(diff0)) {
      ++ np;
      cout << "\n";
      cout << "\ndiff= " << diff;
      cout << "\n";
      for(int ig=0; ig < ndata; ig++) {
        if(flag[ig]) cout << ig << " ";
      }
      cout << "\n";
      for(int ig=0; ig < ndata; ig++) {
        if(flag[ig]) cout << fracs[ig] << " ";
      }
    }
    // Fill the histogram
    hDiff->Fill(diff);
    ++ntrials; ++irun;

  } while( std::next_permutation(flag.begin(), flag.end()) );

  hDiff->Draw();
  // Result
  cout << "np/ntrials = " << double(np) << " / " << double(ntrials) << " = " << double(np)/double(ntrials) << endl;

  return 0;
}
