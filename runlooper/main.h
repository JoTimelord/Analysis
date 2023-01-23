// RAPIDO
#include "arbol.h"
#include "hepcli.h"
#include "looper.h"
#include "cutflow.h"
// ROOT
#include "TString.h"
#include "Math/VectorUtil.h"
// NanoCORE
#include "Nano.h"
#include "Config.h"
#include "tqdm.h"
// Custom
#include "cuts.h"

// Writes the necessary arbol branches such as mass of VVH system
void initializeArbol(Arbol& arbol);
// Holds the "unflattened" variables such as 'LorentzVector' of electrons
void initializeCutflow(Cutflow& cutflow_);

