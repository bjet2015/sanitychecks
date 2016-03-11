#include <functional>
#include <map>

//#define dict map<TString,float>

//list-absed dict - 66 s (use this one)
//my stupid map - 66 s
//std::map - 110 s...
class dict
{
  vector<TString> keys;
  vector<float> values;
public:
  void insert(TString key, float value)
  {
    auto p = std::find(keys.begin(), keys.end(), key);
    if (p == keys.end()) {
      keys.push_back(key);
      values.push_back(value);
    }
    else values[p-keys.begin()]=value;
  }
  //this method can insert in principle, but doesn't do that
  //to make sure values are not spoiled during reading
  bool brokenGet = false;
  float &operator[](const TString key)
  {
    auto p = std::find(keys.begin(), keys.end(), key);
    if (p == keys.end()) {
      cout<<"Key \""<<key<<"\" not found!"<<endl;
      brokenGet = true;
    }
    
    return values[p-keys.begin()];

  }

};

class mystupidmap
{
  const static int maxN = 1<<16;
  TString keys[maxN];
  float values[maxN];
public:
  float &operator[](const TString key)
  {
    unsigned h = key.Hash() % maxN;//hashstring(key);
    keys[h] = key;
    return values[h];
  }
};


bool fillportion = false;

void Fill(TFile *f, vector<TString> varsNeeded, const std::function<void(dict &)> & func)
{
  for (auto x:varsNeeded) cout<<"\""<<x<<"\",";
  cout<<endl;
	TTreeReader reader("nt",f);
	vector<TTreeReaderValue<float> *> values (varsNeeded.size());
	for (unsigned i=0;i<varsNeeded.size();i++)
		values[i] = new TTreeReaderValue<float>(reader,varsNeeded[i]);

	dict v;
	int nev = reader.GetEntries(true);
	int onep = nev/100;
	int evCounter = 0;
	TTimeStamp t0;

	while (reader.Next()) {
	  //if (fillportion && evCounter>2*onep) break;
		evCounter++;
		if (evCounter%onep==0) {
			std::cout << std::fixed; TTimeStamp t1; 
			cout<<" \r"<<evCounter/onep<<"%   "<<" total time "<<(int)round((t1-t0)*nev/(evCounter+.1))<<" s "<<flush;
		}
float a = 0;

		for (unsigned i=0;i<varsNeeded.size();i++)
		  v.insert(varsNeeded[i], *(*(values[i])));//[(TString)varsNeeded[i]] = *(*(values[i]));
		func(v);
    if (v.brokenGet) {
      cout<<"You tried to access key which is not in the list of required branches"<<endl;
      break;
    }
	}
	cout<<endl;
}


vector<TString> concat(vector<TString> a, vector<TString> b)
{
	auto v = vector<TString>();
	v.insert(v.end(),a.begin(), a.end());
	v.insert(v.end(),b.begin(), b.end());
	return v;
}
 
