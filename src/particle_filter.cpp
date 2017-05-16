/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    // This line creates a normal (Gaussian) distribution for x
    
    num_particles = 75; //TODO: Decide how this number is selected.
    default_random_engine gen;
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);
    
    particles.resize(num_particles); //reserve memory for these
    weights.resize(num_particles);
    for (int i = 0; i < num_particles; ++i) {
        particles[i].id = i;
        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen);
        particles[i].weight = 1/num_particles;

    }
    is_initialized = true;
}


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    default_random_engine gen;

    normal_distribution<double> dist_x(0, std_pos[0]);
    normal_distribution<double> dist_y(0, std_pos[1]);
    normal_distribution<double> dist_theta(0, std_pos[2]);
    
    for (int i = 0; i < num_particles; ++i) {
        //Now update the values
        if(fabs(yaw_rate) <= 0.001) {
            particles[i].x += velocity*cos(particles[i].theta)*delta_t;
            particles[i].y += velocity*sin(particles[i].theta)*delta_t;
        } else {
            double theta_next = particles[i].theta + yaw_rate * delta_t;
            particles[i].x += velocity/yaw_rate * (sin(theta_next) - sin(particles[i].theta));
            particles[i].y += velocity/yaw_rate * (cos(particles[i].theta) - cos(theta_next));
            particles[i].theta += yaw_rate*delta_t;
        }
        particles[i].x += dist_x(gen);
        particles[i].y += dist_y(gen);
        particles[i].theta += dist_theta(gen);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

    double total_weight = 0;
    double std_x = std_landmark[0];
    double std_y = std_landmark[1];
    double std_x2 = 2*pow(std_x,2);
    double std_y2 = 2*pow(std_y,2);
    double term_std_xy = (1/(2*M_PI*std_x*std_y)); //calculate outside for loop for speed
    
    for (int i = 0; i < num_particles; ++i) {
        double x = particles[i].x;
        double y = particles[i].y;
        double theta = particles[i].theta;
        
        double prob = 1.0; //Initialized to 1 so we can multipy to itself.
        //Find map landmarks that are with-in sensor_range so they can be detected
        //only for those landmrks look which observations are matching
        Map::single_landmark_s map_lmark;
        for (Map::single_landmark_s landmark : map_landmarks.landmark_list) {
            if(dist(x, y, landmark.x_f, landmark.y_f) > sensor_range) {
                continue;
            }
            double min_dist = numeric_limits<double>::max(); //set to max
            LandmarkObs matched_obs;
            matched_obs.id = -1;
            for (LandmarkObs observation : observations) {
                //transform observation to map coordinates from vehicle coordinates
                LandmarkObs obs;
                obs.x = x + observation.x * cos(theta) - observation.y * sin(theta);
                obs.y = y + observation.x * sin(theta) + observation.y * cos(theta);
            
                //find out closet observation to this landmark.
                double cur_dist = dist(obs.x, obs.y, landmark.x_f, landmark.y_f);
                if (cur_dist <= min_dist) {
                    matched_obs = obs; //store the matched observation
                    matched_obs.id = landmark.id_i;
                    min_dist = cur_dist;
                }
            }
            //calculate Multivariate-Gaussian Probability
            if(matched_obs.id != -1) {
                double d_x2 = pow((matched_obs.x - landmark.x_f), 2)/std_x2; //took stdx_2 calc outside of for loop to speed up.
                double d_y2 = pow((matched_obs.y - landmark.y_f), 2)/std_y2;
                
                double temp_prob = term_std_xy * exp(-1*(d_x2 + d_y2));
                prob *= temp_prob;
                //cout << "found temp_prob " << d_x2 <<"," << d_y2 << endl;
                //cout << "found temp_prob " << temp_prob << endl;
            }
        }
        //cout << "found temp_prob " << prob << endl;
        //store the unnormalized probability as weight for each particle
        if(prob != 1.0) {
            particles[i].weight = prob;
            total_weight += prob;
            //cout << "found new prob " << total_weight << endl;
        } else {
            particles[i].weight = 0; //we didn't find any landmark
        }
    }
    //Now normalize all weights
    for (int i = 0; i < num_particles; ++i) {
        particles[i].weight = particles[i].weight / total_weight;
        weights[i] = particles[i].weight;
        //cout << "weights[i] = " << weights[i] << endl;
    }
     
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    std::vector<Particle> particles_curr(particles);
    std::default_random_engine gen;
    std::discrete_distribution<std::size_t> dis(weights.begin(), weights.end());
    
    for (int i = 0; i < num_particles; ++i) {
        particles[i] = particles_curr[dis(gen)];
    }
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
