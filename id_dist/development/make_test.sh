#!/bin/sh

TEST_DIR=../test
ARGS="-B 1000000 this.file.contains"

rm -f $TEST_DIR/*.f

cp idd_r_test.f $TEST_DIR
cp idz_r_test.f $TEST_DIR
cp idd_a_test.f $TEST_DIR
cp idz_a_test.f $TEST_DIR

grep $ARGS idd_frm.f > $TEST_DIR/idd_frm_test.f
grep $ARGS idd_house.f > $TEST_DIR/idd_house_test.f
grep $ARGS idd_id2svd.f > $TEST_DIR/idd_id2svd_test.f
grep $ARGS idd_id.f > $TEST_DIR/idd_id_test.f
grep $ARGS iddp_aid.f > $TEST_DIR/iddp_aid_test.f
grep $ARGS iddp_asvd.f > $TEST_DIR/iddp_asvd_test.f
grep $ARGS iddp_rid.f > $TEST_DIR/iddp_rid_test.f
grep $ARGS iddp_rsvd.f > $TEST_DIR/iddp_rsvd_test.f
grep $ARGS idd_qrpiv.f > $TEST_DIR/idd_qrpiv_test.f
grep $ARGS iddr_aid.f > $TEST_DIR/iddr_aid_test.f
grep $ARGS iddr_asvd.f > $TEST_DIR/iddr_asvd_test.f
grep $ARGS iddr_rid.f > $TEST_DIR/iddr_rid_test.f
grep $ARGS iddr_rsvd.f > $TEST_DIR/iddr_rsvd_test.f
grep $ARGS idd_sfft.f > $TEST_DIR/idd_sfft_test.f
grep $ARGS idd_snorm.f > $TEST_DIR/idd_snorm_test.f
grep $ARGS idd_svd.f > $TEST_DIR/idd_svd_test.f
grep $ARGS id_rand.f > $TEST_DIR/id_rand_test.f
grep $ARGS id_rtrans.f > $TEST_DIR/id_rtrans_test.f
grep $ARGS idz_frm.f > $TEST_DIR/idz_frm_test.f
grep $ARGS idz_house.f > $TEST_DIR/idz_house_test.f
grep $ARGS idz_id2svd.f > $TEST_DIR/idz_id2svd_test.f
grep $ARGS idz_id.f > $TEST_DIR/idz_id_test.f
grep $ARGS idzp_aid.f > $TEST_DIR/idzp_aid_test.f
grep $ARGS idzp_asvd.f > $TEST_DIR/idzp_asvd_test.f
grep $ARGS idzp_rid.f > $TEST_DIR/idzp_rid_test.f
grep $ARGS idzp_rsvd.f > $TEST_DIR/idzp_rsvd_test.f
grep $ARGS idz_qrpiv.f > $TEST_DIR/idz_qrpiv_test.f
grep $ARGS idzr_aid.f > $TEST_DIR/idzr_aid_test.f
grep $ARGS idzr_asvd.f > $TEST_DIR/idzr_asvd_test.f
grep $ARGS idzr_rid.f > $TEST_DIR/idzr_rid_test.f
grep $ARGS idzr_rsvd.f > $TEST_DIR/idzr_rsvd_test.f
grep $ARGS idz_sfft.f > $TEST_DIR/idz_sfft_test.f
grep $ARGS idz_snorm.f > $TEST_DIR/idz_snorm_test.f
grep $ARGS idz_svd.f > $TEST_DIR/idz_svd_test.f
