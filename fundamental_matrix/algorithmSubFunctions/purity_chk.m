function [ pur ] = purity_chk( VS, s)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

pur =  max(hist(s(VS),max(s))) / length(VS);

