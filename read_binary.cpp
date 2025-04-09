#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <filesystem>
#include <string>
#include <map>
#include <limits>
#include <numbers>
#include <array>
#include <cstdint>
#include "read_binary.hpp"


namespace LSB 
{
    const double Fi    = (99310.0/2147483647.0 * 18.8/4) * (std::numbers::pi/(3600 * 180.00));  // угл.с/имп * pi/(3600.0 * 180.00) = рад/имп
    const double A     = -0.00000003065;
    //const double Tsist =  0.000001; // c
    const double Tdelay=  0.000001; // c
    const double Tsi   =  1.0/48.0;  // мкс
          double T_lg  =  200.0/65535; //0.004;
          double bias_lg = -105;
          double Tadc  = 0.026;
          double bias_adc = -36.00;
          double Ta    =  0.026;
          double bias_a = -36.00;
          double Tsb   =  0.026;
          double bias_sb = -36.00;
}

unsigned int counter = 0;

#pragma pack(1)


Data::Data(Pack & pack): Npack(pack.Npack), Tsi(pack.Tsi * LSB::Tsi), time{0}, Tsist(pack.Tsist), ski(pack.ski), mi(pack.mi),
                         P{0,0,0}, U{0,0,0}, I{0,0,0,0,0,0}, T_lgX{0,0,0}, T_lgY{0,0,0}, T_lgZ{0,0,0}, Tadc{0}, Ta{0,0,0}, Tsb{0}
{
    for (auto i = 0; i < 3; i++)
    {
        A[i] = pack.A[i];
    }

    for (auto i = 0; i < 3; i++)
    {
        Fi[i] = pack.Fi[i];
    }

    for (auto i = 0; i < 3; i++)
    {
        Tdelay[i] = pack.Tdelay[i] * LSB::Tdelay;
    }
}


void Data::count_V(Data & data_old)
{
    for (auto i = 0; i < 3; i++)
    {
        if ((A[i]-data_old.A[i]) > std::numeric_limits<int>::max())
            V.at(i) = ((A[i] - data_old.A[i]) - pow(2,32)) * LSB::A;
        else if ((A[i]-data_old.A[i]) < std::numeric_limits<int>::min())
            V[i] = ((A[i] - data_old.A[i]) + pow(2,32)) * LSB::A;
        else
            V[i] = (A[i] - data_old.A[i]) * LSB::A;
    }
}


void Data::count_dFi(Data & data_old)
{
    for (auto i = 0; i < 3; i++)
    {
        if ((Fi[i]-data_old.Fi[i]) > std::numeric_limits<int>::max())
            dFi.at(i) = ((Fi[i] - data_old.Fi[i]) - pow(2,32)) * LSB::A;
        else if ((Fi[i]-data_old.Fi[i]) < std::numeric_limits<int>::min())
            dFi[i] = ((Fi[i] - data_old.Fi[i]) + pow(2,32)) * LSB::Fi;
        else
            dFi[i] = (Fi[i] - data_old.Fi[i]) * LSB::Fi;
    }
}


void Data::count_W()
{
    for (auto i = 0; i < 3; i++)
    {
        W[i] = V[i] * 1000000 / Tsi;
    }
}

void Data::count_theta()
{
    for (auto i = 0; i < 3; i++)
    {
        Theta[i] = dFi[i] * 1000000 / Tsi;
    }
}


void Data::process_mi()
{
    switch(Npack % 32)
    {
        case 0:
            T_lgX[0] = mi * LSB::T_lg + LSB::bias_lg;
            break;
        case 1:
            T_lgX[1] = mi * LSB::T_lg + LSB::bias_lg;
            break;
        case 2:
            T_lgX[2] = mi * LSB::T_lg + LSB::bias_lg;
            break;
        case 3:
            T_lgY[0] = mi * LSB::T_lg + LSB::bias_lg;
            break;
        case 4:
            T_lgY[1] = mi * LSB::T_lg + LSB::bias_lg;
            break;
        case 5:  
            T_lgY[2] = mi * LSB::T_lg + LSB::bias_lg;
            break;
        case 6:  
            T_lgZ[0] = mi * LSB::T_lg + LSB::bias_lg;
            break;
        case 7: 
            T_lgZ[1] = mi * LSB::T_lg + LSB::bias_lg;
            break;
        case 8: 
            T_lgZ[2] = mi * LSB::T_lg + LSB::bias_lg;
            break;
        case 12:  
            Tadc  = mi * LSB::Tadc + LSB::bias_adc;
            break;
        case 15: 
            Ta[0] = mi * LSB::Ta + LSB::bias_a;
            break;
        case 16: 
            Ta[1] = mi * LSB::Ta + LSB::bias_a;
            break;
        case 17: 
            Ta[2] = mi * LSB::Ta + LSB::bias_a;
            break;
        case 18: 
            Tsb   = mi * LSB::Tsb + LSB::bias_sb;
            break;
        case 19: 
            version = mi;
            break;
        case 20: 
            P[0]  = mi;
            break;
        case 21: 
            P[1]  = mi;
            break;
        case 22:
            P[2]  = mi;
            break;
        case 23: 
            U[0]  = mi;
            break;
        case 24: 
            U[1]  = mi;
            break;
        case 25: 
            U[2]  = mi;
            break;
        case 26: 
            I[0]  = mi;
            break;
        case 27: 
            I[1]  = mi;
            break;
        case 28: 
            I[2]  = mi;
            break;
        case 29: 
            I[3]  = mi;
            break;
        case 30: 
            I[4]  = mi;
            break;
        case 31: 
            I[5]  = mi;
            break;
    } 
}


void Data::process_mi(Data & data)
{
    switch(Npack % 32)
    {
        case 0:
            copy_mi(data);
            T_lgX[0] = mi * LSB::T_lg + LSB::bias_lg;
            break;
        case 1: 
            copy_mi(data);
            T_lgX[1] = mi * LSB::T_lg + LSB::bias_lg;
            break;
        case 2:  
            copy_mi(data);
            T_lgX[2] = mi * LSB::T_lg + LSB::bias_lg;
            break;
        case 3:  
            copy_mi(data);
            T_lgY[0] = mi * LSB::T_lg + LSB::bias_lg;
            break;
        case 4:  
            copy_mi(data);
            T_lgY[1] = mi * LSB::T_lg + LSB::bias_lg;
            break;
        case 5:  
            copy_mi(data);
            T_lgY[2] = mi * LSB::T_lg + LSB::bias_lg;
            break;
        case 6:  
            copy_mi(data);
            T_lgZ[0]= mi * LSB::T_lg + LSB::bias_lg;
            break;
        case 7: 
            copy_mi(data); 
            T_lgZ[1] = mi * LSB::T_lg + LSB::bias_lg;
            break;
        case 8:  
            copy_mi(data);
            T_lgZ[2] = mi * LSB::T_lg + LSB::bias_lg;
            break;
        case 12:  
            copy_mi(data);
            Tadc  = mi * LSB::Tadc + LSB::bias_adc;
            break;
        case 15: 
            copy_mi(data);
            Ta[0] = mi * LSB::Ta + LSB::bias_a;
            break;
        case 16: 
            copy_mi(data);
            Ta[1] = static_cast<short>(mi) * LSB::Ta + LSB::bias_a;
            break;
        case 17: 
            copy_mi(data);
            Ta[2] = mi * LSB::Ta + LSB::bias_a;
            break;
        case 18: 
            copy_mi(data);
            Tsb = mi * LSB::Tsb + LSB::bias_sb;
            break;
        case 19: 
            copy_mi(data);
            version = mi;
            break;
        case 20: 
            copy_mi(data);
            P[0]  = mi;
            break;
        case 21: 
            copy_mi(data);
            P[1]  = mi;
            break;
        case 22:
            copy_mi(data); 
            P[2]  = mi;
            break;
        case 23: 
            copy_mi(data);
            U[0]  = mi;
            break;
        case 24: 
            copy_mi(data);
            U[1]  = mi;
            break;
        case 25: 
            copy_mi(data);
            U[2]  = mi;
            break;
        case 26: 
            copy_mi(data);
            I[0]  = mi;
            break;
        case 27: 
            copy_mi(data);
            I[1]  = mi;
            break;
        case 28: 
            copy_mi(data);
            I[2]  = mi;
            break;
        case 29: 
            copy_mi(data);
            I[3]  = mi;
            break;
        case 30: 
            copy_mi(data);
            I[4]  = mi;
            break;
        case 31: 
            copy_mi(data);
            I[5]  = mi;
            break;
        default:
            copy_mi(data);
    } 
}


void Data::copy_mi(Data & data)
{
    for (auto i = 0; i < 3; i++)
    {
        P[i]  = data.P[i];
        U[i]  = data.U[i];
        Ta[i] = data.Ta[i];
    }
    for (auto i = 0; i < 6; i++)
    {
        I[i]  = data.I[i];
    }
    Tadc      = data.Tadc;
    Tsb       = data.Tsb;
    T_lgX[0]  = data.T_lgX[0];
    T_lgX[1]  = data.T_lgX[1];
    T_lgX[2]  = data.T_lgX[2];
    T_lgY[0]  = data.T_lgY[0];
    T_lgY[1]  = data.T_lgY[1];
    T_lgY[2]  = data.T_lgY[2];
    T_lgZ[0]  = data.T_lgZ[0];
    T_lgZ[1]  = data.T_lgZ[1];
    T_lgZ[2]  = data.T_lgZ[2];
}


void DataSum::add_data(const Data & data)
{
    T      += data.Tsi;
    Npacks += 1;
    Tsist   = data.Tsist;
    Tsi    += data.Tsi;
    time    = data.time;
    Npack   = data.Npack; 
    
    for (auto i = 0; i < 3; i++)
    {
        Tdelay[i] += data.Tdelay[i];
        V[i]     += data.V[i];
        W[i]     += data.W[i];
        P[i]     += data.P[i];
        dFi[i] += data.dFi[i];
        Ta[i]    += data.Ta[i];
        Theta[i] += data.Theta[i];
    }

    for (auto i = 0; i < 6; i++)
    {
        I[i] = I[i] + data.I[i];
    } 

    T_lgX[0]+= data.T_lgX[0];
    T_lgX[1]+= data.T_lgX[1];
    T_lgX[2]+= data.T_lgX[2];
    T_lgY[0]+= data.T_lgY[0];
    T_lgY[1]+= data.T_lgY[1];
    T_lgY[2]+= data.T_lgY[2];
    T_lgZ[0]+= data.T_lgZ[0];
    T_lgZ[1]+= data.T_lgZ[1];
    T_lgZ[2]+= data.T_lgZ[2];
    Tadc    += data.Tadc;
    Tsb     += data.Tsb;
}


void DataSum_m::add_data(const Data & data, const Constants & c)
{
    T      += data.Tsi;
    Npacks += 1;
    Tsist   = data.Tsist;
    Tsi    += data.Tsi;
    time    = data.time;
    Npack   = data.Npack; 
    
    for (auto i = 0; i < 3; i++)
    {
        Tdelay[i]+= data.Tdelay[i];
        V[i]     += data.V[i];
        W[i]     += data.W[i];
        P[i]     += data.P[i];
        dFi[i] += data.dFi[i];
        Ta[i]    += data.Ta[i];
        Theta[i] += data.Theta[i];
    }

    for (auto i = 0; i < 6; i++)
    {
        I[i] = I[i] + data.I[i];
    } 

    T_lgX[0]+= data.T_lgX[0];
    T_lgX[1]+= data.T_lgX[1];
    T_lgX[2]+= data.T_lgX[2];
    T_lgY[0]+= data.T_lgY[0];
    T_lgY[1]+= data.T_lgY[1];
    T_lgY[2]+= data.T_lgY[2];
    T_lgZ[0]+= data.T_lgZ[0];
    T_lgZ[1]+= data.T_lgZ[1];
    T_lgZ[2]+= data.T_lgZ[2];
    Tadc    += data.Tadc;
    Tsb     += data.Tsb;

    std::array<double,3> data_Theta_corr = correct_angles(c, data);
    std::array<double,3> data_dFi_corr;
    std::array<double,3> data_V_corr = correct_V(c, data);
    std::array<double,3> data_W_corr;

    for (auto i = 0; i < 3; i++)
    {
        data_dFi_corr[i] = (data.Tsi/1000000) * data_Theta_corr[i];
        data_W_corr[i] = data_V_corr[i] * 1000000 / data.Tsi;
    }

    for (auto i = 0; i < 3; i++)
    {
        dFi_corr[i] += data_dFi_corr[i];
        Theta_corr[i] += data_Theta_corr[i];
        V_corr[i] += data_V_corr[i];
        W_corr[i] += data_W_corr[i];
    }
}


unsigned short Crc16(unsigned char* data, unsigned short length) 
{
    unsigned short crc = 0xFFFF;
    for (size_t i = 0; i < length; i++) 
    {
        crc ^= data[i];
        for (int j = 0; j < 8; j++) {
            if (crc & 0x0001) 
            { 
                crc >>= 1;      
                crc ^= 0xA001; 
            }
            else 
            {
                crc >>= 1;      
            }
        }
    }
    return crc;
}


void find_header(std::fstream & fin)
{
        //проверяем, неполный ли пакет
    unsigned short x;

    for (auto pos = 47; pos != 0; pos--) //позиция отн начала плохого пакета 
    {
        fin.seekp(-1, std::ios::cur);
        fin.read(reinterpret_cast<char *>(&x), 2);
        fin.seekp(-2, std::ios::cur);
        if (x == 50625)   // C5C1 = 50625 header
        {
            //std::cout << "After pack " << counter << " incomplete pack" << std::endl;
            return;
        }
    }
        //ищем заголовок после плосле плохого пакета
    fin.seekp(48, std::ios::cur); //пришли в конец плохого пакета +1

    while (fin) //пока файл не кончится
    {
        fin.read(reinterpret_cast<char *>(&x), 2);
        if (x == 50625)
            return;
    }
    //больше нет пакетов!!!!!! файл кончился!!!!!
    return;
}
//void find_header(std::fstream & fin)


bool is_good(std::fstream & fin)
{
    char           data[46];
    unsigned short header;
    unsigned short sum;

    fin.read(reinterpret_cast<char*>(&header), 2); 
    fin.seekp(-2, std::ios::cur);
    fin.read(data, 46);

    if (fin.eof())  // проверка на конец файла
    {
        return 0;
    }
    fin.read(reinterpret_cast<char*>(&sum), 2); 
    
    if (header != 50625)
    {
        //вывод в файл с ошибками!!!!!!!!!!!!!!!
        if ((fin.tellg())!=48)   //если нули в начале файла
            std::cout << "After pack " << counter << " header error" << std::endl;
        find_header(fin);
        return 0;
    }
    if (sum != Crc16(reinterpret_cast<unsigned char *>(data), 46))
    {
        unsigned short next;
        fin.read(reinterpret_cast<char *>(&next), 2);
        fin.seekp(-2, std::ios::cur);

        if (next == 50625)
        {
            //целый пакет но битый!!!!!!
            std::cout << "After pack " << counter << " control sum error" << std::endl;
            return 0;
        }
        
        else
        {
            //после пакета не идет сразу следующий
            find_header(fin); //перемещаем курсор к следующему заголовку
            return 0;
        }
    }
    return 1;
}
//bool is_good(std::fstream & fin)


void check_conditions(unsigned short c)
{
    //int nbit = 0; // - номер бита
    if (c & 1) // бит D0
    {
        // равен 1
    }
    c = c >> 1;

    if (c & 1) // бит D1
    {}
    c = c >> 1;

    if (c & 1) // бит D2
    {}
    c = c >> 1;

    if (c & 1) // бит D3
    {}
    c = c >> 1;

    if (c & 1) // бит D4
    {}
    c = c >> 1;

    if (c & 1) // бит D5
    {}
    c = c >> 1;

    if (c & 1) // бит D6
    {}
    c = c >> 1;

    if (c & 1) // бит D7
    {}
    c = c >> 1;

    if (c & 1) // бит D8
    {}
    c = c >> 1;

    if (c & 1) // бит D9
    {}
    c = c >> 1;

    if (c & 1) // бит D10
    {}
    c = c >> 1;

    if (c & 1) // бит D11
    {}
    c = c >> 1;

    if (c & 1) // бит D12
    {}
    c = c >> 1;

    if (c & 1) // бит D13
    {}
    c = c >> 1;

    if (c & 1) // бит D14
    {}

}
//void check_conditions(unsigned short c)


void print_header(std::fstream & fout, std::map<std::string, bool> & output_flags)
{
    if (output_flags["decod"])
    {
        fout << std::setw(10) << "Time[s]";
        fout << std::setw(7)  << "Npack";
        fout << std::setw(10) << "Tsi[mks]";

        fout << std::setw(16) << "dFi_x[rad]"
             << std::setw(16) << "dFi_y[rad]"
             << std::setw(16) << "dFi_z[rad]";
    
        fout << std::setw(14) << "v_x[mps]"
             << std::setw(14) << "v_y[mps]"
             << std::setw(14) << "v_z[mps]";
    
        fout << std::setw(10)  << "Ta1[gC]"
             << std::setw(10)  << "Ta2[gC]"
             << std::setw(10)  << "Ta3[gC]"; 

        fout << std::setw(10)  << "TG11[gC]"
             << std::setw(10)  << "TG12[gC]"
             << std::setw(10)  << "TG13[gC]"
             << std::setw(10)  << "TG21[gC]"
             << std::setw(10)  << "TG22[gC]"
             << std::setw(10)  << "TG23[gC]"
             << std::setw(10)  << "TG31[gC]"
             << std::setw(10)  << "TG32[gC]"
             << std::setw(10)  << "TG33[gC]";
    
        fout << std::setw(10)  << "Tadc[gC]";
        fout << std::setw(10)  << "Tsb[gC]";
        fout << std::setw(6)   << "ski";
        fout << std::endl;
    }
    else if (output_flags["decod_corr"])
    {
        fout << std::setw(10) << "Time[s]";
        fout << std::setw(7)  << "Npack";
        fout << std::setw(10) << "Tsi[mks]";

        fout << std::setw(16) << "dFi_x_corr[rad]"
             << std::setw(16) << "dFi_y_corr[rad]"
             << std::setw(16) << "dFi_z_corr[rad]";
    
        fout << std::setw(14) << "v_x_corr[mps]"
             << std::setw(14) << "v_y_corr[mps]"
             << std::setw(14) << "v_z_corr[mps]";
    
        fout << std::setw(10)  << "Ta1[gC]"
             << std::setw(10)  << "Ta2[gC]"
             << std::setw(10)  << "Ta3[gC]"; 

        fout << std::setw(10)  << "TG11[gC]"
             << std::setw(10)  << "TG12[gC]"
             << std::setw(10)  << "TG13[gC]"
             << std::setw(10)  << "TG21[gC]"
             << std::setw(10)  << "TG22[gC]"
             << std::setw(10)  << "TG23[gC]"
             << std::setw(10)  << "TG31[gC]"
             << std::setw(10)  << "TG32[gC]"
             << std::setw(10)  << "TG33[gC]";
    
        fout << std::setw(10)  << "Tadc[gC]";
        fout << std::setw(10)  << "Tsb[gC]";
        fout << std::setw(6)   << "ski";
        fout << std::endl;
    }
    else if (output_flags["bins"])
    {
        fout << std::setw(16) << "dFi_x_corr[rad]"
             << std::setw(16) << "dFi_y_corr[rad]"
             << std::setw(16) << "dFi_z_corr[rad]";
    
        fout << std::setw(14) << "v_x_corr[mps]"
             << std::setw(14) << "v_y_corr[mps]"
             << std::setw(14) << "v_z_corr[mps]";
        fout << std::endl;
    }
    else
    {
        if (output_flags["Time"])
            fout << std::setw(10) << "Time[s]";
        if (output_flags["Tsist"])
            fout << std::setw(10) << "Tsist[s]";
        if (output_flags["Tsi"])
            fout << std::setw(10) << "Tsi[mks]";
        if (output_flags["Npack"])
            fout << std::setw(7)  << "Npack";
        if (output_flags["A"])
        {
            fout << std::setw(14) << "AK1" 
                 << std::setw(14) << "AK2"
                 << std::setw(14) << "AK3";
        }
        if (output_flags["M"])
        {
            fout << std::setw(16) << "M0"
                 << std::setw(16) << "M1"
                 << std::setw(16) << "M2"
                 << std::setw(16) << "M3";
        }
        if (output_flags["L"])
        {   
            fout << std::setw(16) << "L0"
                 << std::setw(16) << "L1"
                 << std::setw(16) << "L2"
                 << std::setw(16) << "L3";
        }
        if(output_flags["dFi"])
        {
            fout << std::setw(16) << "dFi_x[rad]"
                 << std::setw(16) << "dFi_y[rad]"
                 << std::setw(16) << "dFi_z[rad]";
        }
        if ((output_flags["dFi_corr"])&&(output_flags["Model"]))
        {
            fout << std::setw(16) << std::right << "dFi_x_corr"
                 << std::setw(16) << "dFi_y_corr"
                 << std::setw(16) << "dFi_z_corr";
        }
        if(output_flags["Theta"])
        {
            fout << std::setw(16) << "Theta_x[rps]"
                 << std::setw(16) << "Theta_y[rps]"
                 << std::setw(16) << "Theta_z[rps]";
        }
        if ((output_flags["Theta_corr"])&&(output_flags["Model"]))
        {
            fout << std::setw(16) << "Theta_x_corr"
                 << std::setw(16) << "Theta_y_corr"
                 << std::setw(16) << "Theta_z_corr";
        }
        if (output_flags["Omega"])
        {
            fout << std::setw(17) << "omega_x[gph]"
                 << std::setw(17) << "omega_y[gph]"
                 << std::setw(17) << "omega_z[gph]";
        }
        if ((output_flags["Omega_corr"])&&(output_flags["Model"]))
        {
            fout << std::setw(17) << "omega_x_corr"
                 << std::setw(17) << "omega_y_corr"
                 << std::setw(17) << "omega_z_corr";
        }
        if (output_flags["V"])
        {
            fout << std::setw(14) << "v_x_[mps]"
                 << std::setw(14) << "v_y_[mps]"
                 << std::setw(14) << "v_z_[mps]";
        }
        if ((output_flags["V_corr"])&&(output_flags["Model"]))
        {
            fout << std::setw(14) << "v_x_corr"
                 << std::setw(14) << "v_y_corr"
                 << std::setw(14) << "v_z_corr";
        }
        if (output_flags["W"])
        {    
            fout << std::setw(14) << "w_x[mpss]"
                 << std::setw(14) << "w_y[mpss]"
                 << std::setw(14) << "w_z[mpss]";
        }
        if ((output_flags["W_corr"])&&(output_flags["Model"]))
        {    
            fout << std::setw(14) << "w_x_corr"
                 << std::setw(14) << "w_y_corr"
                 << std::setw(14) << "w_z_corr";
        }
        if (output_flags["Ta"])
        {
            fout << std::setw(10) << "Ta1_[gC]"
                 << std::setw(10) << "Ta2_[gC]"
                 << std::setw(10) << "Ta3_[gC]"; 
        }
        if (output_flags["T_lg"])
        {
            fout << std::setw(10) << "TlgX1[gC]"
                 << std::setw(10) << "TlgX2[gC]"
                 << std::setw(10) << "TlgX3[gC]"
                 << std::setw(10) << "TlgY1[gC]"
                 << std::setw(10) << "TlgY2[gC]"
                 << std::setw(10) << "TlgY3[gC]"
                 << std::setw(10) << "TlgZ1[gC]"
                 << std::setw(10) << "TlgZ2[gC]"
                 << std::setw(10) << "TlgZ3[gC]";
        }
        if (output_flags["T_lg0"])
        {
            fout << std::setw(10) << "TlgX[gC]"
                 << std::setw(10) << "TlgY[gC]"
                 << std::setw(10) << "TlgZ[gC]";
        }
        if (output_flags["Tadc"])
            fout << std::setw(10)  << "Tadc[gC]";
        if (output_flags["Tsb"])
            fout << std::setw(10) << "Tsb[gC]";
        if (output_flags["ski"])
            fout << std::setw(6) << "ski";
        if (output_flags["P"])
        {
            fout << std::setw(9) << "P_1"
                 << std::setw(9) << "P_2"
                 << std::setw(9) << "P_3";
        }
        if (output_flags["U"])
        {
            fout << std::setw(9) << "U_1"
                 << std::setw(9) << "U_2"
                 << std::setw(9) << "U_3";
        }
        if (output_flags["I"])
        {
            fout << std::setw(9) << "I_1"
                 << std::setw(9) << "I_2"
                 << std::setw(9) << "I_3"
                 << std::setw(9) << "I_4"
                 << std::setw(9) << "I_5"
                 << std::setw(9) << "I_6";
        }
        fout << std::endl;
    }
}
//void print_header(std::fstream & fout)

void output_corr_dFi(Data & data, std::fstream & fout, const Constants c)
{
    std::array<double,3> Theta = correct_angles(c, data);
    std::array<double,3> dFi;
    for (auto i = 0; i < 3; i++)
    {
        dFi[i] = Theta[i] * data.Tsi / 1000000;
    }
    fout << std::setw(16) << std::setprecision(12) << std::fixed << dFi[0]
         << std::setw(16) << dFi[1]
         << std::setw(16) << dFi[2];
}


void output_corr_Theta(Data & data, std::fstream & fout, const Constants c)
{
    std::array<double,3> Theta = correct_angles(c, data);
    fout << std::setw(16) << std::setprecision(12) << std::fixed << Theta[0]
         << std::setw(16) << Theta[1]
             << std::setw(16) << Theta[2];
}


void output_corr_Omega(Data & data, std::fstream & fout, const Constants c)
{
    std::array<double,3> Theta = correct_angles(c, data);
    fout << std::setw(17) << std::setprecision(8) << std::fixed << (Theta[0]*180/std::numbers::pi) * 3600
             << std::setw(17) << (Theta[1]*180/std::numbers::pi) * 3600
             << std::setw(17) << (Theta[2]*180/std::numbers::pi) * 3600;
}


void output_corr_V(Data & data, std::fstream & fout, const Constants c)
{
    std::array<double,3> V = correct_V(c, data);
    fout << std::setw(14) << std::setprecision(8) << std::fixed << V[0]
            << std::setw(14) << V[1]
            << std::setw(14) << V[2];
}


void output_corr_W(Data & data, std::fstream & fout, const Constants c)
{
    std::array<double,3> V = correct_V(c, data);
    std::array<double,3> W;
    for (auto i = 0; i < 3; i++)
        W[i] = V[i] * 1000000 / data.Tsi;
    fout << std::setw(14) << std::setprecision(8) << std::fixed << W[0]
            << std::setw(14) << W[1]
            << std::setw(14) << W[2];
}


void output_data(Data & data, std::fstream & fout, std::map<std::string, bool> & output_flags, const Constants & c)
{
    if (output_flags["decod"])
    {
        fout << std::setw(10) << std::setprecision(4)  << std::fixed << data.time;
        fout << std::setw(7)  << data.Npack;
    
        fout << std::setw(10) << std::setprecision(4)  << std::fixed << data.Tsi;
    
        fout << std::setw(16) << std::setprecision(12) << std::fixed << data.dFi[0]
             << std::setw(16) << data.dFi[1]
             << std::setw(16) << data.dFi[2];
    
        fout << std::setw(14) << std::setprecision(8) << std::fixed << data.V[0]
             << std::setw(14) << data.V[1]
             << std::setw(14) << data.V[2];
    
        fout << std::setw(10)  << std::setprecision(2) << std::fixed << data.Ta[0]
             << std::setw(10)  << data.Ta[1]
             << std::setw(10)  << data.Ta[2]; 
    
        fout << std::setw(10)  << std::setprecision(2) << std::fixed << data.T_lgX[0]
             << std::setw(10)  << data.T_lgX[0]
             << std::setw(10)  << data.T_lgX[0]
             << std::setw(10)  << data.T_lgY[0]
             << std::setw(10)  << data.T_lgY[0]
             << std::setw(10)  << data.T_lgY[0]
             << std::setw(10)  << data.T_lgZ[0]
             << std::setw(10)  << data.T_lgZ[0]
             << std::setw(10)  << data.T_lgZ[0];

        fout << std::setw(10)  << std::setprecision(2) << std::fixed << data.Tadc;
        fout << std::setw(10) << (data.T_lgX[0] + data.T_lgY[0] + data.T_lgZ[0])/3;
        fout << std::setw(6) << data.ski;
        fout << std::endl;
    }
    else if (output_flags["decod_corr"])
    {
        fout << std::setw(10) << std::setprecision(4)  << std::fixed << data.time;
        fout << std::setw(7)  << data.Npack;
    
        fout << std::setw(10) << std::setprecision(4)  << std::fixed << data.Tsi;
    
        output_corr_dFi(data, fout, c);
    
        output_corr_V(data, fout, c);
    
        fout << std::setw(10)  << std::setprecision(2) << std::fixed << data.Ta[0]
             << std::setw(10)  << data.Ta[1]
             << std::setw(10)  << data.Ta[2]; 
    
        fout << std::setw(10)  << std::setprecision(2) << std::fixed << data.T_lgX[0]
             << std::setw(10)  << data.T_lgX[0]
             << std::setw(10)  << data.T_lgX[0]
             << std::setw(10)  << data.T_lgY[0]
             << std::setw(10)  << data.T_lgY[0]
             << std::setw(10)  << data.T_lgY[0]
             << std::setw(10)  << data.T_lgZ[0]
             << std::setw(10)  << data.T_lgZ[0]
             << std::setw(10)  << data.T_lgZ[0];

        fout << std::setw(10)  << std::setprecision(2) << std::fixed << data.Tadc;
        fout << std::setw(10) << (data.T_lgX[0] + data.T_lgY[0] + data.T_lgZ[0])/3;
        fout << std::setw(6) << data.ski;
        fout << std::endl;
    }
    else if (output_flags["bins"])
    {
        output_corr_dFi(data, fout, c);
        output_corr_V(data, fout, c);
        fout << std::endl;
    }
    else
    {
        if (output_flags["Time"])
            fout << std::setw(10) << std::setprecision(4)  << std::fixed << data.time;
        if (output_flags["Tsist"])
            fout << std::setw(10) << std::setprecision(4)  << std::fixed << data.Tsist / 1000000.0; 
        if (output_flags["Tsi"])
            fout << std::setw(10) << std::setprecision(4)  << std::fixed << data.Tsi;
        if (output_flags["Npack"])
            fout << std::setw(7)  << data.Npack;
        if (output_flags["dFi"])
        {
            fout << std::setw(16) << std::setprecision(12) << std::fixed << data.dFi[0]
                 << std::setw(16) << data.dFi[1]
                 << std::setw(16) << data.dFi[2];
        }
        if ((output_flags["dFi_corr"])&&(output_flags["Model"]))
            output_corr_dFi(data, fout, c);
        if (output_flags["Theta"])
        {
            fout << std::setw(16) << std::setprecision(12) << std::fixed << data.Theta[0]
                 << std::setw(16) << data.Theta[1]
                 << std::setw(16) << data.Theta[2];
        }
        if ((output_flags["Theta_corr"])&&(output_flags["Model"]))
            output_corr_Theta(data, fout, c);
        if (output_flags["Omega"])
        {
            fout << std::setw(17) << std::setprecision(8) << std::fixed << (data.dFi[0]*180/std::numbers::pi) * 3600 / (data.Tsi/1000000)
                 << std::setw(17) << (data.dFi[1]*180/std::numbers::pi) * 3600 / (data.Tsi/1000000)
                 << std::setw(17) << (data.dFi[2]*180/std::numbers::pi) * 3600 / (data.Tsi/1000000);
        }
        if ((output_flags["Omega_corr"])&&(output_flags["Model"]))
            output_corr_Omega(data, fout, c);
        if (output_flags["V"])
        {
            fout << std::setw(14) << std::setprecision(8) << std::fixed << data.V[0]
                 << std::setw(14) << data.V[1]
                 << std::setw(14) << data.V[2];
        }
        if ((output_flags["V_corr"])&&(output_flags["Model"]))
            output_corr_V(data, fout, c);
        if (output_flags["W"])
        {    
            fout << std::setw(14) << std::setprecision(8) << std::fixed << data.W[0]
                 << std::setw(14) << data.W[1]
                 << std::setw(14) << data.W[2];
        }
        if ((output_flags["W_corr"])&&(output_flags["Model"]))
            output_corr_W(data, fout, c);
        if (output_flags["Ta"])
        {
            fout << std::setw(10)  << std::setprecision(2) << std::fixed << data.Ta[0]
                 << std::setw(10)  << data.Ta[1]
                 << std::setw(10)  << data.Ta[2]; 
        }
        if (output_flags["T_lg"])
        {
            fout << std::setw(10)  << std::setprecision(2) << std::fixed << data.T_lgX[0]
                 << std::setw(10)  << data.T_lgX[1]
                 << std::setw(10)  << data.T_lgX[2]
                 << std::setw(10)  << data.T_lgY[0]
                 << std::setw(10)  << data.T_lgY[1]
                 << std::setw(10)  << data.T_lgY[2]
                 << std::setw(10)  << data.T_lgZ[0]
                 << std::setw(10)  << data.T_lgZ[1]
                 << std::setw(10)  << data.T_lgZ[2];
        }
        if (output_flags["T_lg0"])
        {
            fout << std::setw(10)  << std::setprecision(2) << std::fixed << data.T_lgX[0]
                 << std::setw(10)  << data.T_lgY[0]
                 << std::setw(10)  << data.T_lgZ[0];
        }
        if (output_flags["Tadc"])
            fout << std::setw(10)  << std::setprecision(2) << std::fixed << data.Tadc;
        if (output_flags["Tsb"])
            fout << std::setw(10) << (data.T_lgX[0] + data.T_lgY[0] + data.T_lgZ[0])/3;
        if (output_flags["ski"])
            fout << std::setw(6) << data.ski;
        if (output_flags["P"])
        {
            fout << std::setw(9)  << std::setprecision(4) << std::fixed << data.P[0]
                 << std::setw(9)  << data.P[1]
                 << std::setw(9)  << data.P[2];
        }
        if (output_flags["U"])
        {
            fout << std::setw(9)  << std::setprecision(4) << std::fixed << data.U[0]
                 << std::setw(9)  << data.U[1]
                 << std::setw(9)  << data.U[2];
        }
        if (output_flags["I"])
        {
            fout << std::setw(9)  << std::setprecision(4) << std::fixed << data.I[0]
                 << std::setw(9)  << data.I[1]
                 << std::setw(9)  << data.I[2]
                 << std::setw(9)  << data.I[3]
                 << std::setw(9)  << data.I[4]
                 << std::setw(9)  << data.I[5];
        }
        if (output_flags["mi"])
        {
            fout << std::setw(9)  << data.mi;
        }
        if (output_flags["mis"])
        {
            fout << std::setw(9)  << static_cast<short>(data.mi);
        }
        fout << std::endl;
    }
}
//void output_data(Data & data, std::fstream & fout)


void count_data(Data & data, Data & data_old)
{
    data.time = data_old.time + (data.Tsi/1000000);
    data.count_V(data_old);
    data.count_W();
    data.count_dFi(data_old);
    data.count_theta();
    data.process_mi(data_old);
}


void read_data(std::fstream & fin, std::fstream & fout, std::pair<std::map<std::string, bool>, std::map<std::string, double>> & output_params, const Constants & c)
{
    std::map<std::string, bool> output_flags{output_params.first};
    std::map<std::string, double> time_params{output_params.second};
    bool begin{(!time_params["T_beg"])};
    bool end{false};

    print_header(fout, output_flags);

    Pack pack;
    Data data_old;
    Data data;

    if (is_good(fin)) // читаем первый номер пакета, чтобы потом все нормально шло
    {
        fin.seekp(-46, std::ios::cur); // к номеру пакета
        fin.read(reinterpret_cast<char *>(&pack), 44);
        fin.seekp(2, std::ios::cur); // конец
        counter = pack.Npack;
    }

    data_old = Data(pack);
    data_old.process_mi();

    while (fin)
    {
        if (is_good(fin))
        {
            fin.seekp(-46, std::ios::cur);
            fin.read(reinterpret_cast<char *>(&pack), 44);

            if (pack.Npack != (counter+1))
            { 
                if (counter+1 == pack.Npack-1)
                {
                    std::cout << "Missing pack " << counter + 1 << std::endl;
                }
                else if (!((pack.Npack - counter) > (std::numeric_limits<unsigned short>::max()/2)))
                    std::cout << "Missing packs " <<  counter + 1 << " - " << pack.Npack - 1 << std::endl;
            }

            data = Data(pack);
            count_data(data, data_old);
            counter = pack.Npack;
            if (!begin)
                begin = (data.time >= time_params["T_beg"]);
            if (begin)
                output_data(data, fout, output_flags, c);
            fin.seekp(2, std::ios::cur);
            data_old = data;

            if (time_params["T_end"])
                end = (data.time >= time_params["T_end"]);
            if (end)
                break;
        }
    }
    std::cout << "File is processed" << std::endl;
}
//void read_data(std::fstream & fin, std::fstream & fout)


void output_average_data(DataSum & datasum, std::fstream & fout, std::map<std::string, bool> & output_flags)  // без модели
{
    if (output_flags["Time"])
        fout << std::setw(10) << std::setprecision(4)  << std::fixed << datasum.time;
    if (output_flags["Tsist"])
        fout << std::setw(10) << std::setprecision(4)  << std::fixed << datasum.Tsist / 1000000;
    if (output_flags["Tsi"])
        fout << std::setw(10) << std::setprecision(4)  << std::fixed << datasum.Tsi / datasum.Npacks;
    if (output_flags["Npack"])
        fout << std::setw(7)  << datasum.Npack;
    if(output_flags["dFi"])
    {
        fout << std::setw(16) << std::setprecision(12) << std::fixed << datasum.dFi[0] / datasum.Npacks
             << std::setw(16) << datasum.dFi[1] / datasum.Npacks
             << std::setw(16) << datasum.dFi[2] / datasum.Npacks;
    }
    if(output_flags["Theta"])
    {
        fout << std::setw(16) << std::setprecision(12) << std::fixed << datasum.Theta[0] / datasum.Npacks
             << std::setw(16) << datasum.Theta[1] / datasum.Npacks
             << std::setw(16) << datasum.Theta[2] / datasum.Npacks;
    }
    if (output_flags["Omega"])                                                
    {
        fout << std::setw(17) << std::setprecision(8) << std::fixed << ((datasum.dFi[0]*180/std::numbers::pi) * 3600 / datasum.Npacks) / ((datasum.Tsi / datasum.Npacks)/1000000)
             << std::setw(17) << ((datasum.dFi[1]*180/std::numbers::pi) * 3600 / datasum.Npacks) / ((datasum.Tsi / datasum.Npacks)/1000000)
             << std::setw(17) << ((datasum.dFi[2]*180/std::numbers::pi) * 3600 / datasum.Npacks) / ((datasum.Tsi / datasum.Npacks)/1000000);
    }
    if (output_flags["V"])
    {
        fout << std::setw(14) << std::setprecision(8) << std::fixed << datasum.V[0] / datasum.Npacks
             << std::setw(14) << datasum.V[1] / datasum.Npacks
             << std::setw(14) << datasum.V[2] / datasum.Npacks;
    }
    if (output_flags["W"])
    {    
        fout << std::setw(14) << std::setprecision(8) << std::fixed << datasum.W[0] / datasum.Npacks
             << std::setw(14) << datasum.W[1] / datasum.Npacks
             << std::setw(14) << datasum.W[2] / datasum.Npacks;
    }
    if (output_flags["Ta"])
    {
        fout << std::setw(10)  << std::setprecision(3) << std::fixed << datasum.Ta[0] / datasum.Npacks
             << std::setw(10)  << datasum.Ta[1] / datasum.Npacks
             << std::setw(10)  << datasum.Ta[2] / datasum.Npacks; 
    }
    if (output_flags["T_lg"])
    {
        fout << std::setw(10)  << std::setprecision(3) << std::fixed << datasum.T_lgX[0] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgX[1] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgX[2] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgY[0] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgY[1] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgY[2] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgZ[0] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgZ[1] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgZ[2] / datasum.Npacks;
    }
    if (output_flags["T_lg0"])
    {
        fout << std::setw(10)  << std::setprecision(3) << std::fixed << datasum.T_lgX[0] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgY[0] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgZ[0] / datasum.Npacks;
    }
    if (output_flags["Tadc"])
        fout << std::setw(10)  << std::setprecision(3) << std::fixed << datasum.Tadc / datasum.Npacks;
    if (output_flags["Tsb"])
        fout << std::setw(10) << (datasum.Tsb)/(datasum.Npacks);
    if (output_flags["P"])
    {
        fout << std::setw(9)  << std::setprecision(4) << std::fixed << datasum.P[0] / datasum.Npacks
             << std::setw(9)  << datasum.P[1] / datasum.Npacks
             << std::setw(9)  << datasum.P[2] / datasum.Npacks;
    }
    if (output_flags["U"])
    {
        fout << std::setw(9)  << std::setprecision(4) << std::fixed << datasum.U[0] / datasum.Npacks
             << std::setw(9)  << datasum.U[1] / datasum.Npacks
             << std::setw(9)  << datasum.U[2] / datasum.Npacks;
    }
    if (output_flags["I"])
    {
        fout << std::setw(9)  << std::setprecision(4) << std::fixed << datasum.I[0] / datasum.Npacks
             << std::setw(9)  << datasum.I[1] / datasum.Npacks
             << std::setw(9)  << datasum.I[2] / datasum.Npacks
             << std::setw(9)  << datasum.I[3] / datasum.Npacks
             << std::setw(9)  << datasum.I[4] / datasum.Npacks
             << std::setw(9)  << datasum.I[5] / datasum.Npacks;
    }
    fout << std::endl;
}
//void output_average_data(DataSum & datasum, std::fstream & fout)


void read_average_data(std::fstream & fin, std::fstream & fout, std::pair<std::map<std::string, bool>, std::map<std::string, double>> & output_params)  // без модели
{
    std::map<std::string, bool> output_flags{output_params.first};
    std::map<std::string, double> time_params{output_params.second};
    double taver{time_params["Taverage"]*1000000};
    bool begin{!time_params["T_beg"]};
    bool end{false};

    Pack pack;
    Data data_old;
    Data data;

    if (is_good(fin)) // читаем первый номер пакета, чтобы потом все нормально шло
    {
        fin.seekp(-46, std::ios::cur); // к номеру пакета
        fin.read(reinterpret_cast<char *>(&pack), 44);
        fin.seekp(2, std::ios::cur); // конец
        counter = pack.Npack;
    }

    data_old = Data(pack);

    if (taver < data_old.Tsi)
    {
        fin.seekp(-48, std::ios::cur);
        const Constants c;
        read_data(fin, fout, output_params, c);
        return;
    }

    print_header(fout, output_flags);

    DataSum datasum{taver};

    while (fin)
    {
        if (is_good(fin))
        {
            fin.seekp(-46, std::ios::cur);
            fin.read(reinterpret_cast<char *>(&pack), 44);

            if (pack.Npack != (counter+1))
            { 
                if (counter+1 == pack.Npack-1)
                {
                    std::cout << "Missing pack " << counter + 1 << std::endl;
                }
                else if (!((pack.Npack - counter) > (std::numeric_limits<unsigned short>::max()/2)))
                    std::cout << "Missing packs " <<  counter + 1 << " - " << pack.Npack - 1 << std::endl;
            }

            data = Data(pack);
            count_data(data, data_old);
            counter = pack.Npack;

            if (!begin)
                begin = ((data.time > time_params["T_beg"])&&(abs(data.time-time_params["T_beg"])>(data.Tsi/2000000)));
            if (begin)
                datasum.add_data(data);

            if (datasum.T >= datasum.Taver)
            {
                output_average_data(datasum, fout, output_flags);
                datasum = DataSum(taver);
            }

            fin.seekp(2, std::ios::cur);
            data_old = data;

            if (time_params["T_end"])
                end = ((datasum.time >= (time_params["T_end"]))||(abs(data.time-time_params["T_end"])<(data.Tsi/2000000)));
            if (end)
                break;
        }
    }
    if (datasum.T != 0)
        output_average_data(datasum, fout, output_flags);

    std::cout << "File is processed" << std::endl;
}
//void read_average_data(std::fstream & fin, std::fstream & fout, double taver)


void output_average_data(DataSum_m & datasum, std::fstream & fout, std::map<std::string, bool> & output_flags)  // с моделью
{
    if (output_flags["Time"])
        fout << std::setw(10) << std::setprecision(4)  << std::fixed << datasum.time;
    if (output_flags["Tsist"])
        fout << std::setw(10) << std::setprecision(4)  << std::fixed << datasum.Tsist / 1000000;
    if (output_flags["Tsi"])
        fout << std::setw(10) << std::setprecision(4)  << std::fixed << datasum.Tsi / datasum.Npacks;
    if (output_flags["Npack"])
        fout << std::setw(7)  << datasum.Npack;
    if(output_flags["dFi"])
    {
        fout << std::setw(16) << std::setprecision(12) << std::fixed << datasum.dFi[0] / datasum.Npacks
             << std::setw(16) << datasum.dFi[1] / datasum.Npacks
             << std::setw(16) << datasum.dFi[2] / datasum.Npacks;
    }
    if (output_flags["dFi_corr"])
    {
        fout << std::setw(16) << std::setprecision(12) << std::fixed << datasum.dFi_corr[0] / datasum.Npacks
             << std::setw(16) << datasum.dFi_corr[1] / datasum.Npacks
             << std::setw(16) << datasum.dFi_corr[2] / datasum.Npacks;
    }
    if(output_flags["Theta"])
    {
        fout << std::setw(16) << std::setprecision(12) << std::fixed << datasum.Theta[0] / datasum.Npacks
             << std::setw(16) << datasum.Theta[1] / datasum.Npacks
             << std::setw(16) << datasum.Theta[2] / datasum.Npacks;
    }
    if(output_flags["Theta_corr"])
    {
        fout << std::setw(16) << std::setprecision(12) << std::fixed << datasum.Theta_corr[0] / datasum.Npacks
             << std::setw(16) << datasum.Theta_corr[1] / datasum.Npacks
             << std::setw(16) << datasum.Theta_corr[2] / datasum.Npacks;
    }
    if (output_flags["Omega"])                                                
    {
        fout << std::setw(17) << std::setprecision(8) << std::fixed << ((datasum.dFi[0]*180/std::numbers::pi) * 3600 / datasum.Npacks) / ((datasum.Tsi / datasum.Npacks)/1000000)
             << std::setw(17) << ((datasum.dFi[1]*180/std::numbers::pi) * 3600 / datasum.Npacks) / ((datasum.Tsi / datasum.Npacks)/1000000)
             << std::setw(17) << ((datasum.dFi[2]*180/std::numbers::pi) * 3600 / datasum.Npacks) / ((datasum.Tsi / datasum.Npacks)/1000000);
    }
    if (output_flags["Omega_corr"])                                                
    {
        fout << std::setw(17) << std::setprecision(8) << std::fixed << ((datasum.dFi_corr[0]*180/std::numbers::pi) * 3600 / datasum.Npacks) / ((datasum.Tsi / datasum.Npacks)/1000000)
             << std::setw(17) << ((datasum.dFi_corr[1]*180/std::numbers::pi) * 3600 / datasum.Npacks) / ((datasum.Tsi / datasum.Npacks)/1000000)
             << std::setw(17) << ((datasum.dFi_corr[2]*180/std::numbers::pi) * 3600 / datasum.Npacks) / ((datasum.Tsi / datasum.Npacks)/1000000);
    }
    if (output_flags["V"])
    {
        fout << std::setw(14) << std::setprecision(8) << std::fixed << datasum.V[0] / datasum.Npacks
             << std::setw(14) << datasum.V[1] / datasum.Npacks
             << std::setw(14) << datasum.V[2] / datasum.Npacks;
    }
    if (output_flags["V_corr"])
    {
        fout << std::setw(14) << std::setprecision(8) << std::fixed << datasum.V_corr[0] / datasum.Npacks
             << std::setw(14) << datasum.V_corr[1] / datasum.Npacks
             << std::setw(14) << datasum.V_corr[2] / datasum.Npacks;
    }
    if (output_flags["W"])
    {    
        fout << std::setw(14) << std::setprecision(8) << std::fixed << datasum.W[0] / datasum.Npacks
             << std::setw(14) << datasum.W[1] / datasum.Npacks
             << std::setw(14) << datasum.W[2] / datasum.Npacks;
    }
    if (output_flags["W_corr"])
    {    
        fout << std::setw(14) << std::setprecision(8) << std::fixed << datasum.W_corr[0] / datasum.Npacks
             << std::setw(14) << datasum.W_corr[1] / datasum.Npacks
             << std::setw(14) << datasum.W_corr[2] / datasum.Npacks;
    }
    if (output_flags["Ta"])
    {
        fout << std::setw(10)  << std::setprecision(3) << std::fixed << datasum.Ta[0] / datasum.Npacks
             << std::setw(10)  << datasum.Ta[1] / datasum.Npacks
             << std::setw(10)  << datasum.Ta[2] / datasum.Npacks; 
    }
    if (output_flags["Tadc"])
        fout << std::setw(10)  << std::setprecision(3) << std::fixed << datasum.Tadc / datasum.Npacks;
    if (output_flags["T_lg"])
    {
        fout << std::setw(10)  << std::setprecision(3) << std::fixed << datasum.T_lgX[0] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgX[1] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgY[0] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgY[1] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgZ[0] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgZ[1] / datasum.Npacks;
    }
    if (output_flags["T_lg0"])
    {
        fout << std::setw(10)  << std::setprecision(3) << std::fixed << datasum.T_lgX[0] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgY[0] / datasum.Npacks
             << std::setw(10)  << datasum.T_lgZ[0] / datasum.Npacks;
    }
    if (output_flags["Tsb"])
        fout << std::setw(10) << (datasum.T_lgX[0] + datasum.T_lgY[0] + datasum.T_lgZ[0])/(3*datasum.Npacks);
    if (output_flags["P"])
    {
        fout << std::setw(9)  << std::setprecision(4) << std::fixed << datasum.P[0] / datasum.Npacks
             << std::setw(9)  << datasum.P[1] / datasum.Npacks
             << std::setw(9)  << datasum.P[2] / datasum.Npacks;
    }
    if (output_flags["U"])
    {
        fout << std::setw(9)  << std::setprecision(4) << std::fixed << datasum.U[0] / datasum.Npacks
             << std::setw(9)  << datasum.U[1] / datasum.Npacks
             << std::setw(9)  << datasum.U[2] / datasum.Npacks;
    }
    if (output_flags["I"])
    {
        fout << std::setw(9)  << std::setprecision(4) << std::fixed << datasum.I[0] / datasum.Npacks
             << std::setw(9)  << datasum.I[1] / datasum.Npacks
             << std::setw(9)  << datasum.I[2] / datasum.Npacks
             << std::setw(9)  << datasum.I[3] / datasum.Npacks
             << std::setw(9)  << datasum.I[4] / datasum.Npacks
             << std::setw(9)  << datasum.I[5] / datasum.Npacks;
    }
    fout << std::endl;
}
//void output_average_data(DataSum & datasum, std::fstream & fout)


void read_average_data(std::fstream & fin, std::fstream & fout, std::pair<std::map<std::string, bool>, std::map<std::string, double>> & output_params, const Constants & c)  // с моделью
{
    std::map<std::string, bool> output_flags{output_params.first};
    std::map<std::string, double> time_params{output_params.second};
    double taver{time_params["Taverage"]*1000000};
    bool begin{!time_params["T_beg"]};
    bool end{false};

    Pack pack;
    Data data_old;
    Data data;

    if (is_good(fin)) // читаем первый номер пакета, чтобы потом все нормально шло
    {
        fin.seekp(-44, std::ios::cur); // к номеру пакета
        fin.read(reinterpret_cast<char *>(&pack), 42);
        fin.seekp(2, std::ios::cur); // конец
        counter = pack.Npack;
    }

    data_old = Data(pack);

    if (taver < data_old.Tsi)
    {
        fin.seekp(-46, std::ios::cur);
        read_data(fin, fout, output_params, c);
        return;
    }

    print_header(fout, output_flags);

    DataSum_m datasum{taver};

    while (fin)
    {
        if (is_good(fin))
        {
            fin.seekp(-44, std::ios::cur);
            fin.read(reinterpret_cast<char *>(&pack), 42);

            if (pack.Npack != (counter+1))
            { 
                if (counter+1 == pack.Npack-1)
                {
                    std::cout << "Missing pack " << counter + 1 << std::endl;
                }
                else if (!((pack.Npack - counter) > (std::numeric_limits<unsigned short>::max()/2)))
                    std::cout << "Missing packs " <<  counter + 1 << " - " << pack.Npack - 1 << std::endl;
            }

            data = Data(pack);
            count_data(data, data_old);
            counter = pack.Npack;

            if (!begin)
                begin = ((data.time > time_params["T_beg"])&&(abs(data.time-time_params["T_beg"])>(data.Tsi/2000000)));
            if (begin)
                datasum.add_data(data, c);

            if (datasum.T >= datasum.Taver)
            {
                output_average_data(datasum, fout, output_flags);
                datasum = DataSum_m(taver);
            }

            fin.seekp(2, std::ios::cur);
            data_old = data;

            if (time_params["T_end"])
                end = ((datasum.time >= (time_params["T_end"]))||(abs(data.time-time_params["T_end"])<(data.Tsi/2000000)));
            if (end)
                break;
        }
    }
    if (datasum.T != 0)
        output_average_data(datasum, fout, output_flags);

    std::cout << "File is processed" << std::endl;
}
//void read_average_data(std::fstream & fin, std::fstream & fout, double taver)

void get_flag(std::string & line, std::map<std::string, bool> & flags)
{
    std::istringstream iss(line);
    bool value;
    std::string name;
    std::getline(iss, name, '=');
    iss >> value;
    flags.emplace(std::make_pair(name, value));
} 

double get_number(std::string & line)
{
    std::istringstream iss(line);
    double value;
    std::string name;
    std::getline(iss, name, '=');
    iss >> value;
    return value;
}

std::pair<std::map<std::string, bool>, std::map<std::string, double>> read_config(std::fstream & config)
{
    double value;
    std::map<std::string, bool> output_flags;
    std::map<std::string, double> time_params;
    std::string line;

    std::getline(config, line); // форматирование
    std::getline(config, line); // ---------  

    std::getline(config, line); // decod
    get_flag(line, output_flags); 

    std::getline(config, line);
    get_flag(line, output_flags); // decod_corr

    std::getline(config, line);
    get_flag(line, output_flags); // bins

    std::getline(config, line); // --------- 
    std::getline(config, line); // параметры модели
    std::getline(config, line); // ---------    

    std::getline(config, line);
    get_flag(line, output_flags); // model

    std::getline(config, line); // --------- 
    std::getline(config, line); // временные параметры
    std::getline(config, line); // ---------    

    std::getline(config, line); // Taverage=
    std::istringstream iss(line);
    std::string name;
    std::getline(iss, name, '=');
    iss >> value;
    time_params.emplace(std::make_pair(name, value));

    std::getline(config, line); // T_beg=
    std::istringstream iss1(line);
    std::string name1;
    std::getline(iss1, name1, '=');
    iss1 >> value;
    time_params.emplace(std::make_pair(name1, value));

    std::getline(config, line); // T_end=
    std::istringstream iss2(line);
    std::string name2;
    std::getline(iss2, name2, '=');
    iss2 >> value;
    time_params.emplace(std::make_pair(name2, value));

    std::getline(config, line); // ----------
    std::getline(config, line); // параметры перевода по температурам
    std::getline(config, line); // ----------

    std::getline(config, line); 
    value = get_number(line);
    if (value!=0) LSB::T_lg = value;

    std::getline(config, line); 
    value = get_number(line);
    if (value!=0) LSB::bias_lg = value;

    std::getline(config, line); 
    value = get_number(line);
    if (value!=0) LSB::Tadc = value;

    std::getline(config, line); 
    value = get_number(line);
    if (value!=0) LSB::bias_adc = value;

    std::getline(config, line); 
    value = get_number(line);
    if (value!=0) LSB::Ta = value;

    std::getline(config, line); 
    value = get_number(line);
    if (value!=0) LSB::bias_a = value;

    std::getline(config, line); 
    value = get_number(line);
    if (value!=0) LSB::Tsb = value;

    std::getline(config, line); 
    value = get_number(line);
    if (value!=0) LSB::bias_sb = value;

    std::getline(config, line); // ----------
    std::getline(config, line); // параметры вывода
    std::getline(config, line); // ----------

    while(std::getline(config, line))
    {
        get_flag(line, output_flags);
    }
    return std::make_pair(output_flags, time_params);
}
//std::pair<std::map<std::string, bool>, double> read_config(std::fstream & config)

std::string str_from_config(std::fstream & config)
{
    std::string line;
    std::getline(config, line);  // pcfd=(file_path=)
    std::istringstream iss(line);
    std::string name;
    std::getline(iss, name, '=');
    std::getline(iss, name, ' ');
    std::getline(config, line); // --------
    return name;
}
