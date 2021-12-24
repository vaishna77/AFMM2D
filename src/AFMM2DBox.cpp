#include "AFMM2DBox.hpp"

FMM2DBox::FMM2DBox () {
  boxNumber		=	-1;
  parentNumber	=	-1;
  for (int l=0; l<4; ++l) {
    childrenNumbers[l]	=	-1;
  }
  for (int l=0; l<8; ++l) {
    neighborNumbers[l]	=	-1;
  }
  for (int l=0; l<16; ++l) {
    innerNumbers[l]		=	-1;
  }
  for (int l=0; l<24; ++l) {
    outerNumbers[l]		=	-1;
  }
}
