# DirectedFlowRun18
This has the picos received from Mike in May during EPD removal

How to run the pico:
git clone https://github.com/prashanthgit/DirectedFlowRun18.git
mkdir DirectedFlowRun18-build
cd DirectedFlowRun18-build
cmake ../DirectedFlowRun18
make StPicoEvent PicoAnalyzer AnalyzePico
#to run:
./AnalyzePico 10 "/Users/sprastar/computing/Run18/27GeV/Au27Au2018picos/27GeV_production/list130.list"
  - 10 is the number of file in the list to analyse
  - list file looks like this
  /Users/sprastar/computing/Run18/27GeV/Au27Au2018picos/27GeV_production/130/st_physics_19130077_raw_4000010.picoDst.root
  /Users/sprastar/computing/Run18/27GeV/Au27Au2018picos/27GeV_production/130/st_physics_19130078_raw_4500005.picoDst.root
