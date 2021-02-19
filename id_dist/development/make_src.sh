#!/bin/sh

SRC_DIR=../src
ARGS="-A 1000000000 this.file.contains"

rm -f $SRC_DIR/*.f

cp dfft.f $SRC_DIR
cp prini.f $SRC_DIR

grep $ARGS idd_frm.f > $SRC_DIR/idd_frm.f
grep $ARGS idd_house.f > $SRC_DIR/idd_house.f
grep $ARGS idd_id2svd.f > $SRC_DIR/idd_id2svd.f
grep $ARGS idd_id.f > $SRC_DIR/idd_id.f
grep $ARGS iddp_aid.f > $SRC_DIR/iddp_aid.f
grep $ARGS iddp_asvd.f > $SRC_DIR/iddp_asvd.f
grep $ARGS iddp_rid.f > $SRC_DIR/iddp_rid.f
grep $ARGS iddp_rsvd.f > $SRC_DIR/iddp_rsvd.f
grep $ARGS idd_qrpiv.f > $SRC_DIR/idd_qrpiv.f
grep $ARGS iddr_aid.f > $SRC_DIR/iddr_aid.f
grep $ARGS iddr_asvd.f > $SRC_DIR/iddr_asvd.f
grep $ARGS iddr_rid.f > $SRC_DIR/iddr_rid.f
grep $ARGS iddr_rsvd.f > $SRC_DIR/iddr_rsvd.f
grep $ARGS idd_sfft.f > $SRC_DIR/idd_sfft.f
grep $ARGS idd_snorm.f > $SRC_DIR/idd_snorm.f
grep $ARGS idd_svd.f > $SRC_DIR/idd_svd.f
grep $ARGS id_rand.f > $SRC_DIR/id_rand.f
grep $ARGS id_rtrans.f > $SRC_DIR/id_rtrans.f
grep $ARGS idz_frm.f > $SRC_DIR/idz_frm.f
grep $ARGS idz_house.f > $SRC_DIR/idz_house.f
grep $ARGS idz_id2svd.f > $SRC_DIR/idz_id2svd.f
grep $ARGS idz_id.f > $SRC_DIR/idz_id.f
grep $ARGS idzp_aid.f > $SRC_DIR/idzp_aid.f
grep $ARGS idzp_asvd.f > $SRC_DIR/idzp_asvd.f
grep $ARGS idzp_rid.f > $SRC_DIR/idzp_rid.f
grep $ARGS idzp_rsvd.f > $SRC_DIR/idzp_rsvd.f
grep $ARGS idz_qrpiv.f > $SRC_DIR/idz_qrpiv.f
grep $ARGS idzr_aid.f > $SRC_DIR/idzr_aid.f
grep $ARGS idzr_asvd.f > $SRC_DIR/idzr_asvd.f
grep $ARGS idzr_rid.f > $SRC_DIR/idzr_rid.f
grep $ARGS idzr_rsvd.f > $SRC_DIR/idzr_rsvd.f
grep $ARGS idz_sfft.f > $SRC_DIR/idz_sfft.f
grep $ARGS idz_snorm.f > $SRC_DIR/idz_snorm.f
grep $ARGS idz_svd.f > $SRC_DIR/idz_svd.f
