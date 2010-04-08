clear all;
clc;

dense_succi=load('height09800_succi.dat');
dense_normal=load('height09800.dat');

%% Plotting data
plot(1:256,dense_succi(:,128),1:256,dense_normal(:,128))

%% Density data
mean(dense_succi(180:220,128))
mean(dense_normal(180:220,128))
mean(dense_succi(50:100,128))
mean(dense_normal(50:100,128))