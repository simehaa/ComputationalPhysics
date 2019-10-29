#include "write_file.h"

void write_new_file(cube positions, vec time_vec, vector<planet> planets_list)
{
  // assuming matrix with the form positions(3, N, j);
  int number_of_planets = planets_list.size();
  int N = time_vec.n_elem;

  for (int j=0; j<number_of_planets; j++) {
    string filename = "./data/planet";
    filename.append(to_string(j));
    filename.append(".txt");

    ofstream ofile;
    ofile.open(filename, ofstream::out | ofstream::trunc);
    // write name of planet
    string title = planets_list[j].name;
    // write mass
    ofile << title << " , Mass in solar masses: " << planets_list[j].mass << endl;
    // write header for columns
    ofile << setw(20) << "x: " << setw(20) << "y: " << setw(20) << "z: ";
    ofile << setw(20) << "t: " << endl;
    // write data
    for (int t=0; t<N; t+=100) {
      ofile << setw(20) << setprecision(10) << positions(0, t, j);
      ofile << setw(20) << setprecision(10) << positions(1, t, j);
      ofile << setw(20) << setprecision(10) << positions(2, t, j);
      ofile << setw(20) << setprecision(10) << time_vec(t) << endl;
    }
    cout << "File \"" << filename << "\" written." << endl;
    ofile.close();
  }
}
