#!/bin/bash
echo "Running Cleanup..."
echo "Cleaning up temp data files..."
rm -f /app/temp_data_files/*
echo "cleaning up downloads..."
rm -rf /app/download/*
echo "cleaning up uploads..."
rm -f /app/uploads/*
echo "done cleaning"
