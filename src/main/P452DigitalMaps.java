package main;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;


public class P452DigitalMaps {

    private DigitalMap _mapDN50;
    private DigitalMap _mapN050;


    public P452DigitalMaps() {

        _mapDN50 = new DigitalMap("DN50.TXT", true);
        _mapN050 = new DigitalMap("N050.TXT", true);

    }

    public double GetDN50(double lon, double lat) { return _mapDN50.GetInterpolatedValue(lon, lat);  }

    public double GetN050(double lon, double lat) { return _mapN050.GetInterpolatedValue(lon, lat);  }


    class DigitalMap {
        double[][] _map;
        int _sizeX, _sizeY;
        private double _spacing;

        public DigitalMap(String path, boolean lastAndFirstColumnEqual) {

            String name = "/maps/" + path;
            try {
                InputStream inputStream = getClass().getResourceAsStream(name);

                List<String> lines = new ArrayList<>();

                BufferedReader br = new BufferedReader(new InputStreamReader(inputStream));
                String line;
                while (null != (line = br.readLine())) {
                    lines.add(line);
                }

                _sizeY = lines.size();

                _sizeX = ParseLine(lines.get(0)).length;

                _map = new double[_sizeY][_sizeX];

                for (int i = 0; i < _sizeY; i++) { /* DO */
                    double[] data = ParseLine(lines.get(i));
                    for (int j = 0; j < _sizeX; j++) {
                        _map[i][j] = data[j];
                    }

                }

            } catch (Exception ex) {
                throw new IllegalArgumentException("Could not load map: '" + name + "'");
            }

            if (lastAndFirstColumnEqual) {
                _spacing = 360.0 / (_map[1].length - 1);
            } else {
                _spacing = 360.0 / (_map[1].length);
            }
        }

        private double[] ParseLine(String line) {
            String[] parts = line.trim().split("\\s+");
            ;

            double[] data = new double[parts.length];

            for (int i = 0; i < parts.length; i++) {
                data[i] = Double.parseDouble(parts[i]);
                //System.out.printf(     "parts.length    =  %f\n"  ,data[i]);
            }
            return data;
        }

        public double GetClosestGridPointValue(double longitude, double latitude) {
            double longitudeOffset = longitude + 180.0;
            double latitudeOffset = 90.0 - latitude;
            int latitudeIndex = (int) (latitudeOffset / _spacing);
            int longitudeIndex = (int) (longitudeOffset / _spacing);

            latitudeIndex %= _sizeY;
            longitudeIndex %= _sizeX;


            double val = _map[latitudeIndex][longitudeIndex];
            return val;
        }

        public double GetInterpolatedValue(double longitude, double latitude) {
            double longitudeOffset = longitude;

            if (longitude < 0.0) {
                longitudeOffset = longitude + 360.0;
            }

            double latitudeOffset = 90.0 - latitude;
            int latitudeIndex = (int) (latitudeOffset / _spacing);
            int longitudeIndex = (int) (longitudeOffset / _spacing);
            //System.out.printf(     "latitudeIndex    =  %d\n"  ,latitudeIndex);
            //System.out.printf(     "longitudeIndex    =  %d\n"  ,longitudeIndex);
            double latitudeFraction = (latitudeOffset / _spacing) - latitudeIndex;
            double longitudeFraction = (longitudeOffset / _spacing) - longitudeIndex;

            double value_ul = _map[latitudeIndex][longitudeIndex];
            double value_ur = _map[latitudeIndex][(longitudeIndex + 1) % _sizeX];
            double value_ll = _map[(latitudeIndex + 1) % _sizeY][longitudeIndex];
            double value_lr = _map[(latitudeIndex + 1) % _sizeY][(longitudeIndex + 1) % _sizeX];

            double interpolatedHeight1 = (longitudeFraction * (value_ur - value_ul)) + value_ul;
            double interpolatedHeight2 = (longitudeFraction * (value_lr - value_ll)) + value_ll;
            double interpolatedHeight3 = latitudeFraction * (interpolatedHeight2 - interpolatedHeight1) + interpolatedHeight1;

            return interpolatedHeight3;
        }
    }
}

