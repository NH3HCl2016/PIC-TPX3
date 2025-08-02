/**
 * @file PostProcess.C
 * @brief A simple ROOT macro to process the TPX3 simulation result
 */

{
    // Load required source files
    gROOT->ProcessLine(".L ../source/MCDataToPixeldata3D.C");
    gROOT->ProcessLine(".L ../source/ComptonPolarimeter.C");
    gROOT->ProcessLine(".L ../source/EnergyHistoProcess.C");
    // gROOT->ProcessLine(".L ../source/SavePeakEvents.C");

    // Specify I/O file names
    std::string inputFileDir = "../../Data/Test";
    std::string outputDir = inputFileDir;

    std::string mcDataInputFileName = inputFileDir + "/Run0_RunData.root";
    std::string mcDataOutputFileName = outputDir + "/Run0_PixelData.root";

    std::string expSpecFile = "../../../../../TPX3 Calibration Code/EnergyHist/EnergySpectrum_exp_1332_new.csv";
    std::string expPixelFile = "../../../../../TPX3 Calibration Code/EnergyHist/PixelNumberHist_exp_1332_new.csv";

    Bool_t checkUntriggered = false;

    Double_t roiLeft = 20;
    Double_t roiRight = 40;

    // Specify spectrum and fit range
	// // Na22 - 1274 keV
    // std::pair<Double_t, Double_t> spectrumRange(600.0, 1500.0);
    // std::pair<Double_t, Double_t> fitRange(1200.0, 1350.0);
	// // Co60 - 1332 keV
    // std::pair<Double_t, Double_t> spectrumRange(600.0, 1500.0);
    // std::pair<Double_t, Double_t> fitRange(1260.0, 1500.0);
	// // Cs137 - 662 keV
    // std::pair<Double_t, Double_t> spectrumRange(40.0, 750.0);
    // std::pair<Double_t, Double_t> fitRange(580.0, 720.0);
	// Am241 - 59.5 keV
    std::pair<Double_t, Double_t> spectrumRange(2.0, 80.0);
    std::pair<Double_t, Double_t> fitRange(20.0, 70.0);
	// // Na22 - 511 keV
    // std::pair<Double_t, Double_t> spectrumRange(200.0, 750.0);
    // std::pair<Double_t, Double_t> fitRange(450.0, 560.0);
	// // Co60 - 1173 keV
    // std::pair<Double_t, Double_t> spectrumRange(600.0, 1500.0);
    // std::pair<Double_t, Double_t> fitRange(1140.0, 1260.0);

    ROOT::Math::XYPoint focalCenter(7, 7.5);
    Double_t focalRadius = 5;
    std::pair<ROOT::Math::XYPoint, Double_t> focalRegion = std::make_pair(focalCenter, focalRadius);

    // Specify geometrical and material parameters
    Double_t sensorThickness = 1;
    Double_t xrfThreshold = 0;
    Double_t radius = 0.055 * 2 * TMath::Sqrt(2);
    Double_t minPts = 0;

    std::string outputFilename = std::to_string((roiLeft + roiRight) / 2);

    // MCDataToPixeldata3D(mcDataInputFileName, mcDataOutputFileName, checkUntriggered);
    // SavePeakEvents(mcDataOutputFileName, outputDir, fitRange);
    EnergyHistoProcess(mcDataInputFileName, outputDir, outputFilename, expSpecFile, expPixelFile, fitRange, spectrumRange);
    // ComptonPolarimeter(mcDataOutputFileName, outputDir + "/", radius, minPts, sensorThickness, spectrumRange, xrfThreshold, focalRegion, fitRange);
}
