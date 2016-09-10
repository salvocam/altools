g++ ./bin/contigAlignmentAnalyzer/src/*.cpp -o ./bin/contigAlignmentAnalyzer/Debug/contigAlignmentAnalyzer -lm -std=c++11
g++ ./bin/genicExtractor/src/*.cpp -o ./bin/genicExtractor/Debug/genicExtractor -lm -std=c++11
g++ ./bin/pileupStatistics/src/*.cpp -o ./bin/pileupStatistics/Debug/pileupStatistics -lm -std=c++11
g++ ./bin/variantAnalyzer/src/*.cpp -o ./bin/variantAnalyzer/Debug/variantAnalyzer -lm -std=c++11
g++ ./bin/coverageAnalyzer/src/*.cpp -o ./bin/coverageAnalyzer/Debug/coverageAnalyzer -lm -std=c++11
g++ ./bin/largeDeletionFinder/src/*.cpp -o ./bin/largeDeletionFinder/Debug/largeDeletionFinder -lm -std=c++11
g++ ./bin/slidingAnalysis/src/*.cpp -o ./bin/slidingAnalysis/Debug/slidingAnalysis -lm -std=c++11
g++ ./bin/readsAligner/src/*.cpp -o ./bin/readsAligner/Debug/readsAligner -lm -std=c++11
gcc ./extra/pileupCorr.c -o pileupCorr

