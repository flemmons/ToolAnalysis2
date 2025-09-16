#include "RecoCluster.h"
#include "ANNIERecoObjectTable.h"
#include <algorithm>

RecoCluster::RecoCluster()
{
  serialise=true;
  ANNIERecoObjectTable::Instance()->NewCluster();
}

RecoCluster::~RecoCluster()
{
  ANNIERecoObjectTable::Instance()->DeleteCluster();
  Reset();
}

void RecoCluster::Reset()
{
  int j=fDigitList.size();
  for(int i=0;i<j;i++){
    delete fDigitList[i];
  }
  fDigitList.clear();
}

static bool CompareTimes(RecoDigit *rd1, RecoDigit *rd2)
{
  return ( rd1->GetCalTime() > rd2->GetCalTime() );
}

void RecoCluster::SortCluster()
{
  sort(fDigitList.begin(), fDigitList.end(), CompareTimes);
}

void RecoCluster::AddDigit(RecoDigit* digit)
{
  fDigitList.push_back(digit);
}

RecoDigit* RecoCluster::GetDigit(int n)
{
  return (RecoDigit*)(fDigitList.at(n));
}

int RecoCluster::GetNDigits()
{
  return fDigitList.size();
}

void RecoCluster::SetClusterMode(int cmode)
{
  fClusterMode = cmode;
}

int RecoCluster::GetClusterMode()
{
  return fClusterMode;
}

void RecoCluster::CalcTime() {
	double sum = 0;
	for (RecoDigit* i_digit : fDigitList) {
		sum += i_digit->GetCalTime();
	}
	clusterTime = sum / GetNDigits();

}

void RecoCluster::CalcCharge() {
    double sum=0;
    for (RecoDigit* i_digit : fDigitList) {
        sum += i_digit->GetCalCharge();
    }
    clusterCharge = sum;
}

void RecoCluster::CalcCB() {
    //calculate unfiltered CB
    double total_Q = 0;
    double total_QSquared = 0;
    for (int i = 0; i < fDigitList.size(); i++) {
        //if(unfilteredDigits->at(i)->GetDigitType()==RecoDigit::PMT8inch){
        if (fDigitList.at(i)->GetFilterStatus()) {
            double tube_charge = fDigitList.at(i)->GetCalCharge();
            total_Q += tube_charge;
            total_QSquared += (tube_charge * tube_charge);
            //}
        }
    }
    //FIXME: Need a method to have the 123 be equal to the number of operating detectors
    double charge_balance = sqrt((total_QSquared) / (total_Q * total_Q) - (1. / 123.));

    clusterCB=charge_balance;
}

//Calculate Angular span
// - AS0 is the maximum angle between any two hit PMTs
// - AS1 is the largest angle between the greatest hit by charge and any of the other hits
// - ASC is a charge-weighted version of AS0
void RecoCluster::CalcAS() {

        double max_angle = 0;
        double angle;
        Position i_position, j_position;
        for (int i = 0; i < fDigitList.size(); i++) {
            i_position = fDigitList.at(i)->GetPosition();
            for (int j = 0; j < fDigitList.size(); j++) {
                if (j == i)continue;
                j_position = fDigitList.at(j)->GetPosition();

                angle = j_position.Angle(i_position);
                if (angle > max_angle)max_angle = angle;

            }
        }
        AS0 = max_angle;



        double max_charge = 0;
        double max_angle2 = 0;
        int max_index = 0;
        Position max_position;
        for (int i = 0; i < fDigitList.size(); i++) {
            if (fDigitList.at(i)->GetCalCharge() > max_charge) {
                max_charge = fDigitList.at(i)->GetCalCharge();
                max_index = i;
                max_position = fDigitList.at(i)->GetPosition();
            }
        }


        for (int i = 0; i < fDigitList.size(); i++) {
            if (i == max_index)continue;
            i_position = fDigitList.at(i)->GetPosition();
            angle = i_position.Angle(max_position);
            if (angle > max_angle2)max_angle2 = angle;
        }
        AS1 = max_angle2;


        double max_chargeangle = 0;
        double chargeangle;
        double iq,jq;
        for (int i = 0; i < fDigitList.size(); i++) {
            i_position = fDigitList.at(i)->GetPosition();
            iq=fDigitList.at(i)->GetCalCharge();
            for(int j = 0; j < fDigitList.size(); j++) {
                if(j==i)continue;
                j_position = fDigitList.at(j)->GetPosition();
                jq = fDigitList.at(j)->GetCalCharge();
                chargeangle=iq*jq*j_position.Angle(i_position);
                if(chargeangle>max_chargeangle)max_chargeangle=chargeangle;
            }
        }

        ASC=max_chargeangle;

}

//Retrieve Angular Span value
double RecoCluster::GetAS(int mode) {
    if(mode==0)return AS0;
    else if(mode==1)return AS1;
    else return 0;
}

void RecoCluster::CalcParameters() {
    this->CalcTime();
    this->CalcCharge();
    this->CalcCB();
    this->CalcAS();

    spacialAverage=this->calcSpacialAverage();

    //Other needed parameters here
}

void RecoCluster::SetParticle(int pID, int pPDG, double eff, double pur) {
    bestParticleID=pID;
    bestParticlePDG=pPDG;
    efficiency=eff;
    purity=pur;
}

int RecoCluster::calcBestParent() {
    std::map<int,double> ParticletoClusterCharge;
    std::vector<int> DigiParents;
    int parentCandidate;
    int maxIndex=0;
    double digitCharge;
    double maxCharge=0;
    int tempPDG = -5;
    for (RecoDigit* i_digit : fDigitList) {
        cout<<"Digit parent list size: "<<i_digit->GetParents().size()<<endl;
        if (i_digit->GetParents().size() == 0) {
            cout<<"Digit found with no parents!!!\n";
        }
        for(int i=0; i<i_digit->GetParents().size(); i++){
        parentCandidate=i_digit->GetParents().at(i);
        cout<<"ClusterParent candidate ID: "<<parentCandidate<<endl;
        digitCharge=i_digit->GetCalCharge();
        if (ParticletoClusterCharge.count(parentCandidate) == 0) {
            ParticletoClusterCharge.emplace(parentCandidate,digitCharge);
        }
        else {
            ParticletoClusterCharge.at(parentCandidate)+=digitCharge;
        }
        }
    }
    for (std::pair<int, double>apair : ParticletoClusterCharge) {
        if (apair.second > maxCharge) {
            maxCharge=apair.second;
            maxIndex=apair.first;
            
        }
    }

    /*for (std::pair<int, int> apair : *fMCParticleIndexMap) {
        if (apair.second == maxIndex) {
            tempPDG=apair.first;
        }
    }*/

    //if(tempPDG==-5) return false;

    cout << "Particle " << maxIndex << " wins with charge " << maxCharge << endl;

    bestParticleID=maxIndex;
    //bestParticlePDG=tempPDG;
    purity=maxCharge/clusterCharge;

    return maxIndex;
}

std::vector<RecoDigit*> RecoCluster::convexHull() {

    double tempX, tempY, tempPhi, tempTheta;
    std::vector<double> effX, effY;
    double sumX=0;
    double sumY=0;
    Position pos;
    //std::vector<RecoDigits> hullDigits;
    std::sort(fDigitList.begin(), fDigitList.end());
    for (RecoDigit* adigit : fDigitList) {
        pos = adigit->GetPosition();
        tempX = sin(pos.GetTheta()) * cos(pos.GetPhi());
        tempY = sin(pos.GetTheta()) * sin(pos.GetPhi());
        effX.push_back(tempX);
        effY.push_back(tempY);
        sumX+=tempX;
        sumY+=tempY;
    }

    TwoDCenter.SetX(sumX/fDigitList.size());
    TwoDCenter.SetY(sumY/fDigitList.size());

    TwoDCenter.SetZ(0); //2D projection center; used for circularity test.

    std::vector<int> lower, upper;
    double diffX, diffY;
    double diffprevX, diffprevY;
    double cross;


    // Lower hull
    for (int i = 0; i < fDigitList.size(); i++) {
        if (lower.size() >= 2) {
            diffprevX = effX[lower[lower.size() - 1]] - effX[lower[lower.size() - 2]];
            diffprevY = effY[lower[lower.size() - 1]] - effY[lower[lower.size() - 2]];
            diffX = effX[i] - effX[lower[lower.size() - 1]];
            diffY = effY[i] - effY[lower[lower.size() - 1]];
            cross = diffprevX * diffY - diffprevY * diffX;
        }
        while (lower.size() >= 2 && cross <= 0) {
            //(lower.size() >= 2 && (lower[lower.size() - 1] - lower[lower.size() - 2]).cross(p - lower[lower.size() - 1]) <= 0) {
            lower.pop_back();
            if (lower.size() >= 2) {
                diffprevX = effX[lower[lower.size() - 1]] - effX[lower[lower.size() - 2]];
                diffprevY = effY[lower[lower.size() - 1]] - effY[lower[lower.size() - 2]];
                diffX = effX[i] - effX[lower[lower.size() - 1]];
                diffY = effY[i] - effY[lower[lower.size() - 1]];
                cross = diffprevX * diffY - diffprevY * diffX;
            }

        }
        lower.push_back(i);
        std::cout<<"lowersize,effX,effY,finalcross: "<<lower.size()<<","<<effX[i] <<","<< effY[i]<<","<<cross<<endl;
    }

    // Upper hull
    for (int i = fDigitList.size() - 1; i >= 0; --i) {
        if (upper.size() >= 2) {
            diffprevX = effX[upper[upper.size() - 1]] - effX[upper[upper.size() - 2]];
            diffprevY = effY[upper[upper.size() - 1]] - effY[upper[upper.size() - 2]];
            diffX = effX[i] - effX[upper[upper.size() - 1]];
            diffY = effY[i] - effY[upper[upper.size() - 1]];
            cross = diffprevX * diffY - diffprevY * diffX;
        }
        while (upper.size() >= 2 && cross <= 0) {
            //(upper[upper.size() - 1] - upper[upper.size() - 2]).cross(points[i] - upper[upper.size() - 1]) <= 0) {
            upper.pop_back();
            if (upper.size() >= 2) {
                diffprevX = effX[upper[upper.size() - 1]] - effX[upper[upper.size() - 2]];
                diffprevY = effY[upper[upper.size() - 1]] - effY[upper[upper.size() - 2]];
                diffX = effX[i] - effX[upper[upper.size() - 1]];
                diffY = effY[i] - effY[upper[upper.size() - 1]];
                cross = diffprevX * diffY - diffprevY * diffX;
            }
        }
        upper.push_back(i);
        std::cout << "uppersize,effX,effY,finalcross: " << upper.size() << "," << effX[i] << "," << effY[i] << "," << cross << endl;
    }

    // Remove the last point of each half because it is repeated at the beginning of the other half
    lower.pop_back();
    upper.pop_back();

    // Concatenate lower and upper hull to get the convex hull
    lower.insert(lower.end(), upper.begin(), upper.end());

    for (int i = 0; i < fDigitList.size(); i++) {
        fDigitList.at(i)->ResetFilter();
    }

    for (int i = 0; i < lower.size(); i++) {
        hullDigits.push_back(fDigitList[lower[i]]);
        fDigitList.at(lower[i])->PassFilter();
        std::cout<<"Adding digit index "<<lower[i]<<endl;
    }

    std::cout<<"found hull of size "<<hullDigits.size()<<" compared to "<<fDigitList.size()<<" total digits and lower size: "<<lower.size()<<endl;

    return hullDigits;
}

double RecoCluster::calcSpacialAverage() {
    double aveqX=0;
    double aveqY=0;
    double aveqZ=0;
    double aveqR;
    double digitQ;
    for (RecoDigit* idigit:fDigitList) {
        digitQ = idigit->GetCalCharge();
        aveqX+= digitQ * idigit->GetPosition().X();
        aveqY+= digitQ * idigit->GetPosition().Y();
        aveqZ+= digitQ * idigit->GetPosition().Z();
    }
    aveqX/=clusterCharge;
    aveqY/=clusterCharge;
    aveqZ/=clusterCharge;
    aveqR=pow(aveqX*aveqX + aveqY*aveqY + aveqZ*aveqZ,0.5);

    return aveqR;
}