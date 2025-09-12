/*
* @Author: Wen Jiaxing
* @Date:   2020-03-11 05:21:04
* @Last Modified by:   NH3HCl
* @Last Modified time: 2022-09-16 10:51:43
*/

#include <iostream>
#include <fstream>

#include <TLatex.h>
#include <TLegend.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <TBranch.h>
#include <TStyle.h>
#include <TMath.h>
#include <TF1.h>
#include <TGraph.h>
#include <sys/stat.h>

#include "../../../../../Basic Code/PixelData/3D/include/PixelData3D.hpp"
// #include "../../../../../Basic Code/PlotTrack/3D/include/PlotTrack3D.hpp"

Double_t BkgFunc(Double_t *x, Double_t *par);
Double_t SigFunc(Double_t *x, Double_t *par);
Double_t FitFunc(Double_t *x, Double_t *par);
std::string DoubleToString(Double_t dbNum, Int_t length);
inline Bool_t FileExists(const std::string &name);
void FillExpData(std::string &specFile, std::string &pixelFile, TH1D *spectrum, TH1D *pixelHist, Bool_t readSpectrum = true);

/**
 * @brief Plot the spectrum and histogram for number of triggered pixels from the given input file and save the data
 * @param inputFilename Name of the input file (*.root)
 * @param outputFilepath Path to create the output file
 * @param outputFilename Name of the output file (just the prefix, the file will be saved as *.pdf)
 * @param expSpecFile Name of the measured spectrum file (full name, *.csv)
 * @param expPixelFile Name of the measured pixel number distribution file (full name, *.csv)
 * @param fitRange Upper and lower bounds for the spectrum fit
 * @param spectrumRange Upper and lower bounds for plotting the spectrum
 * @param pixelCountPeakOnly Whether to count only the pixel number of the events in the full energy peak region (center +/- 3 * sigma)
 * @param checkIntvPos Where to check for the local pixel number distribution, in keV
 * @param checkIntvWidth Width of the energy range for local pixel number distribution checking, in keV
 * @param checkExpFiles Name of the files containing the experimental pixel number distributions corresponding to the check energy ranges (*.csv)
 */
void EnergyHistoProcess(std::string inputFilename, std::string outputFilepath, std::string outputFilename, std::string expSpecFile, std::string expPixelFile, std::pair<Double_t, Double_t> fitRange, std::pair<Double_t, Double_t> spectrumRange, Bool_t pixelCountPeakOnly = true, std::vector<Double_t> checkIntvPos = {}, Double_t checkIntvWidth = 0., std::vector<std::string> checkExpFiles = {}) {
    if (FileExists(inputFilename)) {
        TH1D *spectrumHist = new TH1D("Energy Spectrum", "Energy Spectrum pixel", 200, spectrumRange.first, spectrumRange.second);
        TH1D *pixelNumberHist = new TH1D("Pixel number", "Pixel number", 100, 0, 100);
        
        std::cout << "Reading file " << inputFilename << std::endl;

        // Read input data from the specified branch of the tree contained in the input file
        std::vector<Double_t> pixelSignals, pixelCenterX, pixelCenterY, pixelCenterZ;
        std::vector<Int_t> triggered;
        TFile *inputFile = new TFile(inputFilename.c_str());
        TNtuple *pixelData = static_cast<TNtuple*>(inputFile->Get("PixelData3D"));
        if (!pixelData) {
            std::cout << "Could not read the ntuple PixelData3D, unable to proceed." << std::endl;
            return;
        }
        TBranch *signalBranch = pixelData->FindBranch("Signal");
        if (!signalBranch) {
            std::cout << "Could not find the Signal branch on the ntuple, unable to proceed." << std::endl;
            return;
        }
        signalBranch->SetObject(&pixelSignals);
        TBranch *pixelXBranch = pixelData->FindBranch("PixelCenterX");
        if (!pixelXBranch) {
            std::cout << "Could not find the PixelCenterX branch on the ntuple, unable to proceed." << std::endl;
            return;
        }
        pixelXBranch->SetObject(&pixelCenterX);
        TBranch *pixelYBranch = pixelData->FindBranch("PixelCenterY");
        if (!pixelYBranch) {
            std::cout << "Could not find the PixelCenterY branch on the ntuple, unable to proceed." << std::endl;
            return;
        }
        pixelYBranch->SetObject(&pixelCenterY);
        TBranch *pixelZBranch = pixelData->FindBranch("PixelCenterZ");
        if (!pixelZBranch) {
            std::cout << "Could not find the PixelCenterZ branch on the ntuple, unable to proceed." << std::endl;
            return;
        }
        pixelZBranch->SetObject(&pixelCenterZ);
        TBranch *triggeredBranch = pixelData->FindBranch("Triggered");
        if (!triggeredBranch) {
            std::cout << "Could not find the Triggered branch on the ntuple, unable to proceed." << std::endl;
            return;
        }
        triggeredBranch->SetObject(&triggered);

        Double_t meanPixelNumber = 0;
        UInt_t eventCount = 0;
        
        // Process each event
        Long64_t entries = pixelData->GetEntries();
        std::cout << "Total "<< entries << " events to be processed"<<std::endl;
        std::vector<std::pair<Double_t, Double_t>> pixelDataStore;
        Int_t plotCount = 0;
        for (Long64_t i = 0; i < entries; i++) {
            if (i % 10000 == 0) {
                // Print progress
                std::cout << i << "/" << entries << " events processed" << std::endl;
                printf("\033[%dA", (1));
            }
            pixelData->GetEntry(i);
            Double_t totalCharge = 0;
            for (Int_t ip = 0; ip < pixelSignals.size(); ip++) {
                if (triggered[ip] > 0) {
                    // Calculate total charge of the event
                    totalCharge += pixelSignals[ip];
                }
            }
            if (fitRange.first < totalCharge && fitRange.second > totalCharge) {
                meanPixelNumber += pixelSignals.size();
                eventCount ++;
                // if ((totalCharge >= 1308.520743 && totalCharge <= 1313.520743 && pixelData3DResult.size() >= 58) || (totalCharge >= 1365.272719 && totalCharge <= 1370.272719 && pixelData3DResult.size() <= 20)) {
                //     std::cout << std::endl << "Event #" << pixelData3DResult[0].GetEventNumber() << ", " << totalCharge << std::endl;
                //     // Plot current track
                //     Double_t pixelPitch = 0.055;
                //     std::vector<ROOT::Math::XYZPoint> pixelCenters;
                //     std::vector<ROOT::Math::XYZVector> pixelWidths;
                //     std::vector<Double_t> pixelCharges;
                //     for (UInt_t j = 0; j < pixelSignals.size(); j++) {
                //         pixelCenters.push_back(ROOT::Math::XYZPoint(pixelCenterX[j], pixelCenterY[j], pixelCenterZ[j]));
                //         pixelWidths.push_back(ROOT::Math::XYZVector(pixelPitch, pixelPitch, 0));
                //         pixelCharges.push_back(pixelSignals[j]);
                //     }
                //     PlotTrack3D *myPlotTrack3D = new PlotTrack3D();
                //     myPlotTrack3D->SetPixels(pixelCenters, pixelWidths, pixelCharges);
                //     TCanvas *cSETA3D = myPlotTrack3D->PlotPixels();
                //     std::string printFileName = "../Track/" + std::to_string(plotCount) + "_Track_" + std::to_string(totalCharge) + ".pdf";
                //     cSETA3D->Print(printFileName.c_str());
                //     delete cSETA3D;
                //     delete myPlotTrack3D;
                //     plotCount++;
                // }
            }
            if (!pixelCountPeakOnly) {
                pixelNumberHist->Fill(pixelSignals.size());
            }
            else {
                pixelDataStore.push_back(std::make_pair(totalCharge, pixelSignals.size()));
            }
            spectrumHist->Fill(totalCharge);
        }
        std::string pixelNumberString = "Mean pixel number = " + std::to_string(Int_t(meanPixelNumber / eventCount));
        std::cout << pixelNumberString << std::endl;

        // Read measuerd data
        TH1D *expSpecHist = new TH1D("Energy Spectrum", "Energy Spectrum pixel", 200, spectrumRange.first, spectrumRange.second);
        TH1D *expPixelHist = new TH1D("Pixel number", "Pixel number", 100, 0, 100);
        FillExpData(expSpecFile, expPixelFile, expSpecHist, expPixelHist);

        // Plot spectrum
        TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
        c1->cd();
        // c1->SetLogy();
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(0);
        // gStyle->SetOptFit(1111);
        spectrumHist->SetTitle("Energy Spectrum");
        spectrumHist->GetYaxis()->SetTitle("Count");
        spectrumHist->GetXaxis()->SetTitle("Energy/keV");
        spectrumHist->SetLineWidth(4);
        // spectrumHist->SetMinimum(1e2);
        spectrumHist->Draw("hist");
        expSpecHist->Draw("same");
        Double_t maxContent = spectrumHist->GetBinContent(spectrumHist->GetMaximumBin());
        // Fit spectrum based on the specified model
        TF1 *fitFunc = new TF1("FitFunc", FitFunc, fitRange.first, fitRange.second, 6);
        // Set the initial values and parameter names
        fitFunc->SetParameters(0, 0, 0, (fitRange.first + fitRange.second) / 2, (fitRange.first + fitRange.second) / 20, 1000);
        fitFunc->SetParNames("p0", "p1", "p2", "Mean_value", "Sigma", "Constant");
        spectrumHist->Fit(fitFunc, "R");
        Double_t *fitPar = fitFunc->GetParameters();
        const Double_t *fitParError = fitFunc->GetParErrors();
        Double_t scaleFactor = (Double_t)spectrumHist->Integral(spectrumHist->FindBin(fitPar[3] - 3 * fitPar[4]), spectrumHist->FindBin(fitPar[3] + 3 * fitPar[4])) / expSpecHist->Integral(expSpecHist->FindBin(fitPar[3] - 3 * fitPar[4]), expSpecHist->FindBin(fitPar[3] + 3 * fitPar[4]));
        expSpecHist->Scale(scaleFactor);
        expPixelHist->Scale(scaleFactor);

        // Create histograms for the given energy regions for pixel number checking
        std::vector<TH1D*> checkHistArray;
        if (checkIntvWidth > 0) {
            for (UInt_t i = 0; i < checkIntvPos.size(); i++) {
                std::string histTitle = "Pixel number in energy range: " + std::to_string(fitPar[3] + checkIntvPos[i] * fitPar[4] - checkIntvWidth / 2) + " - " + std::to_string(fitPar[3] + checkIntvPos[i] * fitPar[4] + checkIntvWidth / 2) + " keV";
                checkHistArray.push_back(new TH1D("Pixel number", histTitle.c_str(), 100, 0, 100));
            }
        }

        // Count the number of triggered pixels in the full energy peak region
        if (pixelCountPeakOnly) {
            for (Int_t i = 1; i < pixelDataStore.size(); i++) {
                if (fitPar[3] - 3. * fitPar[4] < pixelDataStore[i].first && fitPar[3] + 3. * fitPar[4] > pixelDataStore[i].first) {
                    pixelNumberHist->Fill(pixelDataStore[i].second);
                }
                // Count the number of triggered pixels in all check regions
                if (checkIntvWidth > 0) {
                    for (UInt_t j = 0; j < checkIntvPos.size(); j++) {
                        if (fitPar[3] + checkIntvPos[j] * fitPar[4] - checkIntvWidth / 2 < pixelDataStore[i].first && fitPar[3] + checkIntvPos[j] * fitPar[4] + checkIntvWidth / 2 > pixelDataStore[i].first) {
                            checkHistArray[j]->Fill(pixelDataStore[i].second);
                        }
                    }
                }
            }
        }

        // Plot fitted signal and background
        TF1 *fitBkg = new TF1("BkgFunc", BkgFunc, fitRange.first, fitRange.second, 3);
        TF1 *fitSig = new TF1("SigFunc", SigFunc, fitRange.first, fitRange.second, 3);
        for (Int_t j = 0; j < 6; j++) {
            fitFunc->SetParameter(j, fitPar[j]);
            if (j < 3) {
                fitBkg->SetParameter(j, fitPar[j]);
            }
            else {
                fitSig->SetParameter(j - 3, fitPar[j]);
            }
        }
        fitBkg->Draw("same");
        fitBkg->SetLineColor(kGreen);
        fitFunc->Draw("same");
        fitFunc->SetLineColor(kRed);
        // fitSig->Draw("same");
        // fitSig->SetLineColor(kBlue);
        
        Double_t peakIntensity = fitBkg->Eval(fitPar[3]) + fitSig->Eval(fitPar[3]);
        std::string gaussExprString = "#it{f(E)}=#it{a}#frac{1}{#sqrt{2#it{#pi}}#it{#sigma}}exp(#frac{#it{(E-#mu)}^{2}}{2#it{#sigma}^{2}})+#it{p_{0}E^{2}+p_{1}E+p_{2}}";
        std::string parameterString = std::string(" #mu=") + DoubleToString(fitPar[3], 2) + "#pm" + DoubleToString(fitParError[3], 2) + " (keV)" + std::string(" #sigma=") + DoubleToString(fitPar[4], 2) + "#pm" + DoubleToString(fitParError[4], 2) + " (keV)";

        // Add legend and labels
        TLegend *legend1 = new TLegend(0.56, 0.71, 0.87, 0.9);
        legend1->AddEntry(spectrumHist, "Simulation", "l");
        legend1->AddEntry(expSpecHist, "Experiment", "l");
        legend1->SetTextSize(0.035);
        legend1->SetMargin(0.1);
        legend1->Draw();

        spectrumHist->SetMaximum(maxContent * 1.5);
        TLatex latex;
        latex.SetTextSize(0.04);
        latex.DrawLatex(spectrumRange.first, maxContent * 1.4, gaussExprString.c_str());
        latex.DrawLatex(spectrumRange.first, maxContent * 1.1, parameterString.c_str());
        latex.DrawLatex(spectrumRange.first, maxContent * 0.9, pixelNumberString.c_str());
        c1->Print(Form("%s/%s_spectrum.pdf", outputFilepath.c_str(), outputFilename.c_str()));
        
        // Plot histogram for number of triggered pixels
        TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
        c2->cd();
        c2->SetLogy();
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(0);
        pixelNumberHist->SetTitle("Pixel Number");
        pixelNumberHist->GetYaxis()->SetTitle("Count");
        pixelNumberHist->GetXaxis()->SetTitle("Pixel Number");
        pixelNumberHist->SetLineWidth(4);
        pixelNumberHist->Draw("hist");
        expPixelHist->Draw("same");

        // Add legend and labels
        TLegend *legend2 = new TLegend(0.56, 0.71, 0.87, 0.9);
        legend2->AddEntry(pixelNumberHist, "Simulation", "l");
        legend2->AddEntry(expPixelHist, "Experiment", "l");
        legend2->SetTextSize(0.035);
        legend2->SetMargin(0.1);
        legend2->Draw();
        
        c2->Print(Form("%s/%s_pixelnumber.pdf", outputFilepath.c_str(), outputFilename.c_str()));
        
        // Plot histograms for check energy ranges
        if (checkIntvWidth > 0) {
            TH1D *checkPixelHist = new TH1D("Pixel number", "Pixel number", 100, 0, 100);
            for (UInt_t i = 0; i < checkIntvPos.size(); i++) {
                // c2->cd();
                // c2->SetLogy();
                // gStyle->SetOptStat(0);
                // gStyle->SetOptFit(0);
                checkHistArray[i]->GetYaxis()->SetTitle("Count");
                checkHistArray[i]->GetXaxis()->SetTitle("Pixel Number");
                checkHistArray[i]->SetLineWidth(4);
                checkHistArray[i]->Draw("hist");
                // Read measured data from the given experimental data file
                checkPixelHist->Clear();
                std::string dummyString = "";
                FillExpData(dummyString, checkExpFiles[i], 0, checkPixelHist, false);
                checkPixelHist->Scale(scaleFactor);
                checkPixelHist->Draw("same");
                std::string currOutFileName = "check_" + std::to_string(fitPar[3] + checkIntvPos[i] * fitPar[4] - checkIntvWidth / 2) + "_" + std::to_string(fitPar[3] + checkIntvPos[i] * fitPar[4] + checkIntvWidth / 2);
                c2->Print(Form("%s/%s_pixelnumber.pdf", outputFilepath.c_str(), currOutFileName.c_str()));
            }
        }

        // Print spectrum info
        std::ofstream specrtrumOut(Form("%s/EnergySpectrum.csv", outputFilepath.c_str()), std::ios::out | std::ios::ate);
        specrtrumOut << "BinCenter(keV)" << ", " << "Counts" << std::endl;
        for (Int_t i = 0; i < spectrumHist->GetNbinsX(); i++) {
            specrtrumOut << spectrumHist->GetBinCenter(i) << ", " << spectrumHist->GetBinContent(i) << std::endl;
        }
        specrtrumOut.close();
        // Print info of histogram for number of triggered pixels
        std::ofstream pixlenumberOut(Form("%s/PixelNumberHist.csv", outputFilepath.c_str()), std::ios::out | std::ios::ate);
        pixlenumberOut << "PixelNumber" << ", " << "Counts" << std::endl;
        for (Int_t i = 0; i < pixelNumberHist->GetNbinsX(); i++) {
            pixlenumberOut << pixelNumberHist->GetBinCenter(i) << ", " << pixelNumberHist->GetBinContent(i) << std::endl;
        }
        pixlenumberOut.close();
    }
}

/**
 * @brief Function for the quadratic background component in the spectrum
 * @param x Input x values
 * @param par Parameters
 * @return Corresponding values of the background function
 */
Double_t BkgFunc(Double_t *x, Double_t *par) {
    Double_t xx = x[0];
    Double_t p0 = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    return p0 + p1 * xx + p2 * xx * xx;
}

/**
 * @brief Function for the gaussian peak in the spectrum
 * @param x Input x values
 * @param par Parameters
 * @return Corresponding values of the peak function
 */
Double_t SigFunc(Double_t *x, Double_t *par) {
    Double_t xx = x[0];
    Double_t mu = par[0];
    Double_t sigma = par[1];
    //Double_t binwidth = 9;
    Double_t norm = par[2];
    return norm * TMath::Exp(-(xx - mu) * (xx - mu) / 2 / sigma / sigma) / TMath::Sqrt(2 * TMath::Pi()) / sigma;
}

/**
 * @brief Function used to fit the spectrum (gaussian peak + quadratic background)
 * @param x Input x values
 * @param par Parameters
 * @return Corresponding values of the fit function
 */
Double_t FitFunc(Double_t *x, Double_t *par) {
    return BkgFunc(x, par) + SigFunc(x, &par[3]);
}

/**
 * @brief Check if the given file exists
 * @param name Name of the file to be checked
 * @return Whether the file exists
 */
inline Bool_t FileExists(const std::string &name) {
  struct stat buffer;   
  return (stat(name.c_str(), &buffer) == 0); 
}

/**
 * @brief Convert a Double_t to string with a given number of decimal places
 * @param dbNum The Double_t to be converted
 * @param length Nmuber of decimal places to keep
 * @return Converted string
 */
std::string DoubleToString(Double_t dbNum, Int_t length) {
    // By default the initial string will have 6 decimal places after converting with TMath::Nint
    Int_t defaultlength = 6;
    dbNum = Double_t(TMath::Nint(dbNum * 100)) / 100;
    std::string str = std::to_string(dbNum);
    str=str.substr(0, str.size() - (defaultlength - length));
    return str;
}

/**
 * @brief Read spectrum and pixel number distribution from input file (*.csv) and fill the corresponding histograms. MAKE SURE the range and binning of the histograms are the same as the corresponding data in the input file
 * @param specFile File containing the spectrum
 * @param pixelFile File containing the pixel number distribution
 * @param spectrum Histogram storing the spectrum
 * @param pixelHist Histogram storing the pixel number distribution
 * @param readSpectrum Whether to read the spectrum data
 */
void FillExpData(std::string &specFile, std::string &pixelFile, TH1D *spectrum, TH1D *pixelHist, Bool_t readSpectrum = true) {
    Int_t i = 0;
    // Read spectrum
    if (readSpectrum) {
        std::ifstream specInput(specFile, std::ios::in);
        if (!specInput.is_open()) {
            std::cout << "Error: unable to open experimental spectrum file \"" << specFile << "\" or file not found, unable to plot experimental spectrum" << std::endl;
        }
        else {
            std::string strLine;
            // Skip the header
            std::getline(specInput, strLine);
            while (std::getline(specInput, strLine)) {
                if (strLine.empty()) {
                    continue;
                }
                std::istringstream strIn(strLine);
                std::string newVal;
                std::getline(strIn, newVal, ',');
                std::getline(strIn, newVal);
                spectrum->SetBinContent(i, std::atof(newVal.c_str()));
                i++;
            }
            specInput.close();
        }
    }

    // Read pixel number distribution
    i = 0;
    std::ifstream pixelInput(pixelFile, std::ios::in);
    if (!pixelInput.is_open()) {
        std::cout << "Error: unable to open experimental pixel number distribution file \"" << pixelFile << "\" or file not found, unable to plot experimental spectrum" << std::endl;
    }
    else {
        std::string strLine;
        // Skip the header
        std::getline(pixelInput, strLine);
        while (std::getline(pixelInput, strLine)) {
            if (strLine.empty()) {
                continue;
            }
            std::istringstream strIn(strLine);
            std::string newVal;
            std::getline(strIn, newVal, ',');
            std::getline(strIn, newVal);
            pixelHist->SetBinContent(i, std::atof(newVal.c_str()));
            i++;
        }
    }
}
