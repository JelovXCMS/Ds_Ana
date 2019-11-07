// NonPrompt DataWeightOverFONLL


// for pp

weightfunctiongen="GptSampleWeight*FonllGptWeight*(0.55369+1.22995*log(GBAncestorpt)+-0.49766*log(GBAncestorpt)*log(GBAncestorpt)+0.05655*log(GBAncestorpt)*log(GBAncestorpt)*log(GBAncestorpt))";
weightfunctionreco="DPtSampleWeight*FonllPtWeight*(0.55369+1.22995*log(DgenBAncestorpt)+-0.49766*log(DgenBAncestorpt)*log(DgenBAncestorpt)+0.05655*log(DgenBAncestorpt)*log(DgenBAncestorpt)*log(DgenBAncestorpt))";



// for PbPb
weightfunctiongen="GptSampleWeight*FonllGptWeight*centralityNCollWeight*vertexZWeight*(1.43138+0.37192*log(GBAncestorpt)+-0.71770*log(GBAncestorpt)*log(GBAncestorpt)+0.23269*log(GBAncestorpt)*log(GBAncestorpt)*log(GBAncestorpt)+-0.02181*log(GBAncestorpt)*log(GBAncestorpt)*log(GBAncestorpt)*log(GBAncestorpt))";
weightfunctionreco="DPtSampleWeight*FonllPtWeight*centralityNCollWeight*vertexZWeight*(1.43138+0.37192*log(DgenBAncestorpt)+-0.71770*log(DgenBAncestorpt)*log(DgenBAncestorpt)+0.23269*log(DgenBAncestorpt)*log(DgenBAncestorpt)*log(DgenBAncestorpt)+-0.02181*log(DgenBAncestorpt)*log(DgenBAncestorpt)*log(DgenBAncestorpt)*log(DgenBAncestorpt))";

