#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <assert.h>
#include <fstream>
#include <x86intrin.h>


#include "seal/seal.h"
#include "Examples.cpp"
#include "slots.h"
#include "transform.h"
#include "bootstrapping.h"
#include "karatsuba.h"


using namespace std;
using namespace seal;

vector<uint64_t> sk_vector_generator(int n, int slot_count, int N_n);
vector<vector<int>> A_matrix_generator(int slot_count, int n);
Ciphertext Hom_Lin_Trans(GaloisKeys galois_keys[], GaloisKeys gk2, vector<vector<uint64_t>> A, Ciphertext encrypted_matrix, PublicKey public_key, RelinKeys relin_keys,SEALContext context, SecretKey secret_key, SlotRing basic_slot_ring); //Declare Hom_Lin_Trans function

int main()
{
    
    // set up context
	EncryptionParameters parms(scheme_type::bfv);
	std::shared_ptr<SlotRing> slot_ring = p_257_test_parameters();//p_257_test_parameters();
	SlotRing basic_slot_ring = slot_ring->change_exponent(1); // this part here has to be 1 since they assume plaintext space is actually p^r = p
	parms.set_poly_modulus_degree(slot_ring->N());
	parms.set_coeff_modulus(CoeffModulus::BFVDefault(slot_ring->N()));
	parms.set_plain_modulus(slot_ring->prime());// this part here has to be prime since they assume plaintext space is actually p^r = p
	SEALContext context(parms);
    print_parameters(context);
    cout << endl;
	std::cout << "Q bit count is " << context.key_context_data()->total_coeff_modulus_bit_count() << std::endl;
    
	// set up keys
	KeyGenerator keygen(context);
	SecretKey sk = keygen.secret_key();
	PublicKey pk;
	keygen.create_public_key(pk);
	RelinKeys rk;
	keygen.create_relin_keys(rk);

	Encryptor encryptor(context, pk);
	Evaluator evaluator(context);
	Decryptor decryptor(context, sk);
    /*
	// create an example plaintext
	poly data;
	poly_add(data, basic_slot_ring.from_slot_value({ 1 }, 0), basic_slot_ring.R().scalar_mod); //basic_slot_ring.R().scalar_mod.value() is just p^e. Note data[0] is now 255 cuz /257 since .change_exponent(2) in slots and /257 since change_exponent(1) above.
	poly_add(data, basic_slot_ring.from_slot_value({ 8 }, 1), basic_slot_ring.R().scalar_mod); // .from_slot_value({ 8 }, 1) puts value 8 into slot 1 of data. Use {8,1} for 8+x etc.
	poly_add(data, basic_slot_ring.from_slot_value({ 63 }, 2), basic_slot_ring.R().scalar_mod);
	poly_add(data, basic_slot_ring.from_slot_value({ 17 }, 3), basic_slot_ring.R().scalar_mod);
	poly_add(data, basic_slot_ring.from_slot_value({ 12 }, 4), basic_slot_ring.R().scalar_mod);
    
	Plaintext x_plain{ gsl::span<const uint64_t>(data) };
    poly data2(x_plain.data(), x_plain.data() + basic_slot_ring.N()); //NOTE HAVE TO DO THIS DO CONVERT PLAINTEXT TO ENCODED. CANNOT JUST EQUATE AT VECTORS

    cout << "N: "<<slot_ring->N()<<endl;
    cout << "Extracted value: " << basic_slot_ring.extract_slot_value(data2, 0)[0] << endl;

    //for (int i=0;i<512;i++){
    //    cout<<basic_slot_ring.extract_slot_value(data2, i)[0] << ' ';
    //}

	Ciphertext x_enc;
	encryptor.encrypt(x_plain, x_enc);
	std::cout << "encrypted message" << std::endl;
	std::cout << "noise budget is " << decryptor.invariant_noise_budget(x_enc) << " bits" << std::endl;
    
    SlotRing::Rotation rot = basic_slot_ring.rotate(128);
    GaloisKeys gk;
    keygen.create_galois_keys(std::vector{ std::get<0>(rot.galois_elements()), std::get<1>(rot.galois_elements()) }, gk);
    Ciphertext result_enc;
	rot.apply_ciphertext(x_enc, evaluator, gk, result_enc);
    Plaintext result;
	decryptor.decrypt(result_enc, result);
    poly result_poly(result.data(), result.data() + basic_slot_ring.N());//result.data() is the address of the poly. So not too sure why theyre adding it with N()

    for (int j=0;j<512;j++){
        if (basic_slot_ring.extract_slot_value(result_poly, j)[0] == 1){
            cout<< "Index i: "<<j<<endl;
        }
    }
    */


    
	
    
    
    
    // MSEAL BOOTSTRAPPING
    

    // Encoding/Creating plaintext (Secret key//bfvct)
    size_t slot_count = 128; 
    size_t n = slot_count/2;
    size_t N_n = slot_count / n;
    size_t LWE_dim = slot_count*1; //Set to *1 if doing n < slot_count

    // Generate sk vector. Dimension could be anything.
    vector<uint64_t> pod_matrix = sk_vector_generator(n, slot_count, N_n); //Set FIRST entry to n<N or LWE_dim. If set to n, it will compute (sk||sk||...||sk).





    
    //Do bootstrapping on sk(LWE) vector
    poly data_boot;
    for (int i=0; i<slot_count; i++){
        // create an example plaintext
	    poly_add(data_boot, basic_slot_ring.from_slot_value({ pod_matrix[i] }, i), basic_slot_ring.R().scalar_mod); 
    }
    Plaintext x_plain_boot{ gsl::span<const uint64_t>(data_boot) };

    Ciphertext x_enc; //Encrypt plaintext
    encryptor.encrypt(x_plain_boot, x_enc);

    Bootstrapper bootstrapper(context, slot_ring, p_257_test_parameters_digit_extractor(*slot_ring));
	bootstrapper.initialize();

	std::cout << "initialized bootstrapper" << std::endl;

	BootstrappingKey bk;
	bootstrapper.create_bootstrapping_key(sk, bk);

	std::cout << "created bootstrapping key" << std::endl;

	auto start = std::chrono::steady_clock::now();

	Ciphertext in_coeffs;
	bootstrapper.slots_to_coeffs(x_enc, bk, in_coeffs, MemoryManager::GetPool() DEBUG_PASS(sk));

	auto end = std::chrono::steady_clock::now();
	std::cout << "moved to coefficients in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms and " << kswitch_counter << " key-switch operations" << std::endl;
	kswitch_counter = 0;
	start = end;

	Ciphertext noisy_dec;
	bootstrapper.homomorphic_noisy_decrypt(in_coeffs, bk, noisy_dec, MemoryManager::GetPool() DEBUG_PASS(sk));

	end = std::chrono::steady_clock::now();
	std::cout << "homomorphically decrypted in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms and " << kswitch_counter << " key-switch operations" << std::endl;
	start = end;
	kswitch_counter = 0;
    




    //Creating a random matrix A
    vector<vector<int>> A = A_matrix_generator(slot_count, n);//Set SECOND entry to n<N or LWE_dim.

    // Calc A*s naively
    vector<uint64_t> rslt(slot_count, 0ULL);
    cout << "Multiplication of given two matrices is:\n";
    for (int i = 0; i < A.size(); i++) {
            rslt[i] = 0;
                for (int k = 0; k < A[0].size(); k++) {
                    rslt[i] += A[i][k] * pod_matrix[k];
                }
        }
    for (int i=0; i<slot_count;i++){
        cout << rslt[i]%257 << ' ';
    }
    cout <<endl;

    //Performing hom. A*s
    // Initialising variables
    vector<uint64_t> x_plain_data(basic_slot_ring.N() , 0ULL); 
    Plaintext plain_matrix{ gsl::span<const uint64_t>(x_plain_data) }; //Encode 0 plaintext
    Ciphertext summed_rslt;
    encryptor.encrypt(plain_matrix, summed_rslt);

    for (int k=0; k < LWE_dim/slot_count; k++){
        int B_2nd_dim;
        if (LWE_dim == slot_count){
            B_2nd_dim = n;
        } else{ 
            B_2nd_dim = slot_count;
        }

        vector<vector<uint64_t> > B(slot_count, vector<uint64_t>(B_2nd_dim));
        poly data;
        for (int i=0; i<slot_count; i++){
            // create an example plaintext
	        poly_add(data, basic_slot_ring.from_slot_value({ pod_matrix[i+k*slot_count] }, i), basic_slot_ring.R().scalar_mod); 
            for (int j = 0; j < B_2nd_dim; j++) { 
                B[i][j] = A[i][j+k*B_2nd_dim];
            } 
        }
        Plaintext x_plain{ gsl::span<const uint64_t>(data) };

        Ciphertext encrypted_matrix; //Encrypt plaintext
        encryptor.encrypt(x_plain, encrypted_matrix);

        GaloisKeys gk[slot_count];
        int rt = sqrt(n);
        
        for (int i=0;i<rt;i++){
            SlotRing::Rotation rot = basic_slot_ring.rotate(-(i+1)*rt);
            keygen.create_galois_keys(std::vector{ std::get<0>(rot.galois_elements()), std::get<1>(rot.galois_elements()) }, gk[i]);
        }
        GaloisKeys gk2; //containing rotation by 1.
        SlotRing::Rotation rot = basic_slot_ring.rotate(-1);
        keygen.create_galois_keys(std::vector{ std::get<0>(rot.galois_elements()), std::get<1>(rot.galois_elements()) }, gk2);

        auto start2 = std::chrono::steady_clock::now();
        Ciphertext rslt_enc = Hom_Lin_Trans(gk, gk2, B, encrypted_matrix, pk, rk, context, sk, basic_slot_ring);
        auto end2 = std::chrono::steady_clock::now();
        std::cout << "homomorphically decrypted in " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count() << " ms  " << std::endl;

        evaluator.add_inplace(summed_rslt, rslt_enc);
    }
    
    print_line(__LINE__);
    Plaintext plain_result; //Used for testing, decrypting outputs
    cout << "Decrypt and decode result." << endl;
    decryptor.decrypt(summed_rslt, plain_result);
    cout << "    + Noise budget in encrypted_matrix: " << decryptor.invariant_noise_budget(summed_rslt) << " bits" << endl;
    poly result_poly(plain_result.data(), plain_result.data() + basic_slot_ring.N());//result.data() is the address of the poly. So not too sure why theyre adding it with N()

    for (int j=0;j<slot_count;j++){
        cout<<basic_slot_ring.extract_slot_value(result_poly, j)[0]<< ' ';
    }
    
    
	return 0;
}


vector<vector<int>> A_matrix_generator(int slot_count, int n)
{
    std::random_device rd2;     // Only used once to initialise (seed) engine
    std::mt19937 rng2(rd2());    // Random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<unsigned long long int> uni2(0,257); // Guaranteed unbiased
    vector<vector<int>> A;
    for (int i = 0; i < slot_count; i++) {  
        vector<int> temp; 
        for (int j = 0; j < n; j++) { // Set LWE dim to n if doing normal mult
            temp.push_back(uni2(rng2)); 
        } 
        A.push_back(temp); 
    }
    return A;
}

vector<uint64_t> sk_vector_generator(int n, int slot_count, int N_n)
{
    vector<uint64_t> pod_matrix2(n, 0ULL);// SK
    std::random_device rd;     // Only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // Random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<unsigned long long int> uni(0,127); // Guaranteed unbiased
    for (int i = 0; i < n; i++) { 
        unsigned long long int random_integer = uni(rng);
        pod_matrix2[i] = random_integer;
        }
    if (n < slot_count){
        vector<uint64_t> pod_matrix(slot_count, 0ULL); //(SK||SK||...||SK)
        for (int i = 0; i<N_n;i++){
            for (int j = 0; j<n;j++){
                pod_matrix[j+i*n] = pod_matrix2[j]; //Creates bfvct
            }
        }
        return pod_matrix;
    }
    return pod_matrix2;
}

Ciphertext Hom_Lin_Trans(GaloisKeys galois_keys[], GaloisKeys gk2, vector<vector<uint64_t>> A, Ciphertext encrypted_matrix, PublicKey public_key, RelinKeys relin_keys,SEALContext context, SecretKey secret_key, SlotRing basic_slot_ring)
{
    Encryptor encryptor(context, public_key); // Do we need to declare this globally or use pointers instead?
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    // Initialising res_k
    int n = A[0].size();
    cout << "n size: "<< n <<endl;
    size_t slot_count = 128;
    int rt = sqrt(n);
    vector<uint64_t> res_k(basic_slot_ring.N() , 0ULL); 

	Plaintext res_k_plain{ gsl::span<const uint64_t>(res_k) };
    Ciphertext res_k_enc;
    encryptor.encrypt(res_k_plain, res_k_enc);
    Ciphertext res[rt];
    auto start4 = std::chrono::steady_clock::now();
    for (int i = 0; i<rt;i++){ //might need to adjust index to start at 1
        res[i] = res_k_enc;
    }
    auto end4 = std::chrono::steady_clock::now();
    std::cout << "Assigned vector values in " << std::chrono::duration_cast<std::chrono::milliseconds>(end4 - start4).count() << " ms  " << std::endl;
    
    //Initialising bfvct_rot_is
    Ciphertext bfvct_rot[rt]; // Initializing bfvct_rot_is

    Ciphertext bfvct_rot2[rt];
    auto start5 = std::chrono::steady_clock::now();
    auto end5 = std::chrono::steady_clock::now();
    for (int i = 0; i<rt; i++){
            bfvct_rot[i] = encrypted_matrix;
            if (i==1){
                start5 = std::chrono::steady_clock::now();
            } 
            SlotRing::Rotation rot = basic_slot_ring.rotate(-(i+1)*rt);
            if (i==1){
                end5 = std::chrono::steady_clock::now();
                std::cout << "Assigned rot in " << std::chrono::duration_cast<std::chrono::milliseconds>(end5 - start5).count() << " ms  " << std::endl;
            }
            rot.apply_ciphertext(bfvct_rot[i], evaluator, galois_keys[i], bfvct_rot2[i]);
    }
    
    Ciphertext c;
    int ind_ct;
    int ind_a;
    for (int k = 1; k<rt+1; k++){
        for (int i = 1; i<rt+1; i++){
            poly tmp;
            for (int j = 0; j<slot_count;j++){
                ind_ct = (j-k+1) % slot_count;
                ind_a = (j+i*rt) % n;
                if (j==1 && i==1 && k==1){
                    start5 = std::chrono::steady_clock::now();
                } 
                poly_add(tmp, basic_slot_ring.from_slot_value({ A[ind_ct][ind_a] }, j), basic_slot_ring.R().scalar_mod);
                if (j==1 && i==1 && k==1){
                    end5 = std::chrono::steady_clock::now();
                    std::cout << "Assigned add in " << std::chrono::duration_cast<std::chrono::milliseconds>(end5 - start5).count() << " ms  " << std::endl;
                }
            }
            Plaintext tmp_plain{ gsl::span<const uint64_t>(tmp) };
            evaluator.multiply_plain(bfvct_rot2[i-1], tmp_plain, c);
            evaluator.relinearize_inplace(c, relin_keys); 
            evaluator.add_inplace(res[k-1],c);
        }
    }

    SlotRing::Rotation rot = basic_slot_ring.rotate(-1);
    for (int i = 1; i<rt; i++){
        rot.apply_ciphertext(res[rt-i+1-1], evaluator, gk2, c);
        evaluator.add_inplace(res[rt-i-1],c);
    }

    return res[0];
}

