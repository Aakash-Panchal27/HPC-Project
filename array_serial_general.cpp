#include<bits/stdc++.h>
#include<fstream>
// #include<omp.h>
using namespace std;

const int bits_unit = 8;

struct node {
    int frequency;
    unsigned char in;
    node * left;
    node * right;
    node(int freq, unsigned char i)
    {
        frequency = freq;
        in = i;
        left = NULL;
        right = NULL;
    }
   
};

struct comparator{

    bool operator()(node* one, node* two)
    {
        if(one->frequency > two->frequency)
            return true;
        else 
            return false;
    }
};

void obtain_huffmancodes(node* root, string current, unordered_map<unsigned char,string>& encode)
{
    if(root->left == NULL && root->right == NULL)
    {
        encode[root->in] = current;
        return;
    }
    obtain_huffmancodes(root->left, string(current + "0"), encode);
    obtain_huffmancodes(root->right, string(current + "1"), encode);
}

int main()
{
    string filename = "ttest.raw";
    string compressed_file_name = "test_file_compressed.raw";
    string decompressed_file_name = "testop.raw";
    ifstream testing(filename, ios::in | ios::binary);
    unsigned char a;

    // Create frequency table--------------------------------------------------
    auto start = chrono::high_resolution_clock::now();

    int freq[256] = {0};
    while(testing.good()){
        testing.read((char*)&a, sizeof(a));
        freq[a]++;
    }
    auto end = chrono::high_resolution_clock::now();
    double total_time =  
      chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    total_time *= 1e-9;
    cout << "Creating frequency table: " << endl;
    cout << fixed << setprecision(20) << total_time << endl;
    testing.close();
    // Create frequency table--------------------------------------------------

    // for(auto it:freq)
    //  cout << it.first << " " << it.second << endl;

    // Create huffman tree-----------------------------------------------------
    start = chrono::high_resolution_clock::now();

    priority_queue<node*, vector<node*>, comparator> container;

    for(int i = 0; i < 256; i++)
    {
    	if(freq[i] > 0) {
	        node* temp = new node(freq[i], i);
	        container.push(temp);
    	}
    }

    while(container.size() > 1)
    {
        node* left = container.top();
        container.pop();
        
        node * right = container.top();
        container.pop();
        
        node* curr = new node(left->frequency + right->frequency, '0');
        curr->left = left;
        curr->right = right;
        container.push(curr);
    }

    node* root = container.top();
    container.pop();
    node* root_store = root;

    end = chrono::high_resolution_clock::now();
    total_time =  
      chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    total_time *= 1e-9;
    cout << "Creating Huffman Tree: " << endl;
    cout << fixed << setprecision(20) << total_time << endl;
    // Create huffman tree-------------------------------------------------

    // Creating Hashtable of encodes----------------------------------------
    start = chrono::high_resolution_clock::now();

    string ans = "";
    unordered_map<unsigned char,string> encode;
    obtain_huffmancodes(root, ans, encode);

    end = chrono::high_resolution_clock::now();
    total_time =  
      chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    total_time *= 1e-9;
    cout << "Creating hashtable of encodes: " << endl;
    cout << fixed << setprecision(20) << total_time << endl;
    // Creating Hashtable of encodes----------------------------------------

    // Compress--------------------------------------------------
    start = chrono::high_resolution_clock::now();

    short int idx = 7;
    int total_bits = 0;
    ifstream final_read(filename, ios::in | ios::binary);
    ofstream compressed(compressed_file_name, ios::out | ios::binary);

    unsigned char c = 0;
    while(final_read.good()){
        final_read.read((char*)&a, sizeof(a));
        string str = encode[a];
        int len = str.length();
        total_bits += len;
        for(int i = 0; i < len; i++){
            bool bit = (str[i] - '0');
            c |= (bit << idx);
            idx--;
            if(idx < 0) {
                compressed.write((char*)&c, sizeof(c));
                idx = 7;
                // refresh c
                c &= 0;
            }
        }
    }

    if(idx != 7)
        compressed.write((char*)&c, sizeof(c));
    final_read.close();
    compressed.close();

    end = chrono::high_resolution_clock::now();
    total_time =  
      chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    total_time *= 1e-9;
    cout << "Compression: " << endl;
    cout << fixed << setprecision(20) << total_time << endl;
    // Compress-------------------------------------------------
    
    // Idk why
    //total_bits--;

    // Decompress---------------------------------------------------
    start = chrono::high_resolution_clock::now();

    ifstream compressed_final_read(compressed_file_name, ios::in | ios::binary);
    ofstream decompressed(decompressed_file_name, ios::out | ios::binary);
    root = root_store;
    while(compressed_final_read.good() && total_bits > 0){
        compressed_final_read.read((char*)&a, sizeof(a));
        idx = 7;
        while(idx >= 0 && total_bits > 0) {
            if(a & (1 << idx))
                root = root->right;
            else
                root = root->left;
            if(!root->right && !root->left) {
                unsigned char cur = (root->in);
                decompressed.write((char *)&cur, sizeof(cur));
                root = root_store;
            }
            idx--;
            total_bits--;
        }
    }
    compressed_final_read.close();
    decompressed.close();

    end = chrono::high_resolution_clock::now();
    total_time =  
      chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    total_time *= 1e-9;
    cout << "Decompression " << endl;
    cout << fixed << setprecision(20) << total_time << endl;
    // Decompress-----------------------------------------------------------
    return 0;
}